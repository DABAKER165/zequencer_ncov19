# First stage: JDK 11 with modules
FROM debian:stretch-slim as packager

# source JDK distribution names
# update from https://jdk.java.net/java-se-ri/11
ENV JDK_VERSION="11.0.1"
ENV JDK_URL="https://download.java.net/java/GA/jdk11/13/GPL/openjdk-${JDK_VERSION}_linux-x64_bin.tar.gz"
ENV JDK_HASH="7a6bb980b9c91c478421f865087ad2d69086a0583aeeb9e69204785e8e97dcfd"
ENV JDK_HASH_FILE="${JDK_ARJ_FILE}.sha2"
ENV JDK_ARJ_FILE="openjdk-${JDK_VERSION}.tar.gz"
# target JDK installation names
ENV OPT="/opt"
ENV JKD_DIR_NAME="jdk-${JDK_VERSION}"
ENV JAVA_HOME="${OPT}/${JKD_DIR_NAME}"
ENV JAVA_MINIMAL="${OPT}/java-minimal"

# downlodad JDK to the local file
ADD "$JDK_URL" "$JDK_ARJ_FILE"
# COPY openjdk-11.0.1_linux-x64_bin.tar.gz "$JDK_ARJ_FILE"
# verify downloaded file hashsum
RUN { \
        echo "Verify downloaded JDK file $JDK_ARJ_FILE:" && \
        echo "$JDK_HASH $JDK_ARJ_FILE" > "$JDK_HASH_FILE" && \
        sha256sum -c "$JDK_HASH_FILE" ; \
    }

# extract JDK and add to PATH
RUN { \
        echo "Unpack downloaded JDK to ${JAVA_HOME}/:" && \
        mkdir -p "$OPT" && \
        tar xf "$JDK_ARJ_FILE" -C "$OPT" ; \
    }
ENV PATH="$PATH:$JAVA_HOME/bin"

RUN { \
        java --version ; \
        echo "jlink version:" && \
        jlink --version ; \
    }

# build modules distribution
RUN jlink \
    --verbose \
    --add-modules \
        java.base,java.sql,java.naming,java.desktop,java.management,java.security.jgss,java.instrument \
        # java.naming - javax/naming/NamingException
        # java.desktop - java/beans/PropertyEditorSupport
        # java.management - javax/management/MBeanServer
        # java.security.jgss - org/ietf/jgss/GSSException
        # java.instrument - java/lang/instrument/IllegalClassFormatException
    --compress 2 \
    --strip-debug \
    --no-header-files \
    --no-man-pages \
    --output "$JAVA_MINIMAL"

RUN apt update \
 && apt install -y upx-ucl wget python3 python3-pip
RUN apt-get install -y build-essential
RUN apt-get install -y libncurses5-dev libncursesw5-dev zlib1g-dev libbz2-dev liblzma-dev bcftools
# samtools

ADD https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz /v1.3.tar.gz
RUN cd / && tar -zxf v1.3.tar.gz && rm v1.3.tar.gz && \
cd seqtk-1.3 && make
ADD https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip /snpEff_latest_core.zip
RUN apt-get install -y unzip
RUN cd / && unzip snpEff_latest_core.zip

ADD https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.4.jar /VarScan.v2.4.4.jar

ADD https://downloads.sourceforge.net/project/bbmap/BBMap_38.86.tar.gz?ts=1617285740&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2FBBMap_38.86.tar.gz%2Fdownload%3Fuse_mirror%3Dphoenixnap /BBMap_38.86.tar.gz
RUN cd / && tar -zxf BBMap_38.86.tar.gz && rm BBMap_38.86.tar.gz

ADD https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 /htslib-1.12.tar.bz2
RUN cd / && tar -xjf htslib-1.12.tar.bz2 && rm htslib-1.12.tar.bz2 && \
cd htslib-1.12 && \
./configure --prefix=/htslib && \
make && \
make install



# # Second stage, add only our minimal "JRE" distr and our app
FROM debian:stretch-slim

ENV JAVA_HOME=/opt/java-minimal
ENV PATH="$PATH:$JAVA_HOME/bin"

COPY --from=packager "$JAVA_HOME" "$JAVA_HOME"
COPY --from=packager /snpEff/snpEff.jar /bin/snpEff.jar
COPY --from=packager /snpEff/SnpSift.jar /bin/SnpSift.jar
COPY --from=packager /snpEff/snpEff.config /bin/snpEff.config
COPY --from=packager /VarScan.v2.4.4.jar /bin/VarScan.v2.4.4.jar
COPY --from=packager /seqtk-1.3/seqtk /bin/seqtk
COPY --from=packager /bbmap /bin/bbmap
# COPY --from=packager /bbmap/reformat.sh /bin/reformat.sh
# COPY --from=packager /bbmap/bbmerge.sh /bin/bbmerge.sh
# COPY --from=packager /bbmap/bbmap.sh /bin/bbmap.sh
# COPY --from=packager /bbmap/bbduk.sh /bin/bbduk.sh
# COPY --from=packager /bbmap/callvariants.sh /bin/callvariants.sh

# COPY --from=packager /usr/bin/samtools /bin/samtools
COPY --from=packager /usr/bin/bcftools /bin/bcftools
COPY --from=packager /htslib/bin/tabix /bin/tabix
COPY --from=packager /htslib/bin/bgzip /bin/bgzip


RUN apt update \
 && apt install -y python3 python3-pip samtools \
 # bcftools \
 && python3 -m pip --no-cache-dir install snakemake biopython==1.76 pyfasta pandas pysam \
 && apt remove -y gcc libpython3.5-dev dpkg-dev openssh-client \
 && apt autoremove -y

RUN ln -s /bcftools/bin/bcftools /usr/bin/bcftools && \
# ln -s /samtools/bin/samtools /usr/bin/samtools && \
ln -s /bin/bbmap/reformat.sh /usr/bin/reformat.sh && \
ln -s /bin/bbmap/bbmerge.sh /usr/bin/bbmerge.sh && \
ln -s /bin/bbmap/bbmap.sh /usr/bin/bbmap.sh && \
ln -s /bin/bbmap/bbduk.sh /usr/bin/bbduk.sh && \
ln -s /bin/bbmap/callvariants.sh /usr/bin/callvariants.sh && \
ln -s /bin/seqtk /usr/bin/seqtk && \
ln -s /bin/tabix /usr/bin/tabix && \
ln -s /bin/bgzip /usr/bin/bgzip && \
# ln -s /bin/novoindex /usr/bin/novoindex && \
# ln -s /bin/novoalign /usr/bin/novoalign && \
printf '#!/bin/sh\njava -jar /bin/snpEff.jar "$@"' > /usr/bin/snpEff.sh && \
chmod +x /usr/bin/snpEff.sh && \
ln -s /usr/bin/snpEff.sh /usr/bin/snpEff && \
printf '#!/bin/sh\njava -jar /bin/VarScan.v2.4.4.jar "$@"' > /usr/bin/varscan.sh && \
chmod +x /usr/bin/varscan.sh && \
ln -s /usr/bin/varscan.sh /usr/bin/varscan

COPY ref/ /ref
COPY nCoV-2019/ /nCoV-2019
COPY zequencer.sh /zequencer.sh
COPY trim_strict.py /trim_strict.py
COPY nCoV-2019_adapters.fasta /nCoV-2019_adapters.fasta
COPY amplicon_normalization.py /amplicon_normalization.py
COPY amplicon_normalize_from_bam.py /amplicon_normalize_from_bam.py
# COPY zequencer.smk /zequencer.smk