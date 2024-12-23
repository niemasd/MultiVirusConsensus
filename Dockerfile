# Minimal Docker image for MultiVirusConsensus using Ubuntu base
FROM ubuntu:24.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# install MultiVirusConsensus and dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y autoconf automake bzip2 cmake g++ libboost-all-dev libbz2-dev libcurl4-openssl-dev liblzma-dev make python3 wget zlib1g-dev && \

    # install htslib (needed for ViralConsensus)
    wget -qO- "https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2" | tar -xj && \
    cd htslib-* && autoreconf -i && ./configure && make && make install && cd .. && rm -rf htslib-* && \

    # install samtools
    wget -qO- "https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2" | tar -xj && \
    cd samtools-* && ./configure --without-curses && make && make install && cd .. && rm -rf samtools-* && \

    # install minimap2
    wget -qO- "https://github.com/lh3/minimap2/archive/refs/tags/v2.28.tar.gz" | tar -zx && \
    cd minimap2-* && make && chmod a+x minimap2 && mv minimap2 /usr/local/bin/minimap2 && cd .. && rm -rf minimap2-* && \

    # install ViralConsensus
    wget -qO- "https://github.com/niemasd/ViralConsensus/archive/refs/tags/0.0.6.tar.gz" | tar -zx && \
    cd ViralConsensus-* && make && mv viral_consensus /usr/local/bin/viral_consensus && cd .. && rm -rf ViralConsensus-* && \

    # install Google Sparsehash (needed for BioBloom)
    wget -qO- "https://github.com/sparsehash/sparsehash/archive/refs/tags/sparsehash-2.0.4.tar.gz" | tar -zx && \
    cd sparsehash-* && ./configure && make && make install && cd .. && rm -rf sparsehash-* && \

    # install sdsl-lite (needed for BioBloom)
    wget -qO- "https://github.com/simongog/sdsl-lite/releases/download/v2.1.1/sdsl-lite-2.1.1.tar.gz.offline.install.gz" | tar -zx && \
    cd sdsl-lite-* && ./install.sh /usr/local/ && cd .. && rm -rf sdsl-lite-* && \

    # install BioBloom
    wget -qO- "https://github.com/bcgsc/biobloom/releases/download/2.3.5/biobloomtools-2.3.5.tar.gz" | tar -zx && \
    cd biobloomtools-* && sed -i 's/c++11/c++14/g' configure.ac && ./configure && make && make install && cd .. && rm -rf biobloomtools-* && \

    # install MultiVirusConsensus
    wget -q -O /usr/local/bin/MultiVirusConsensus.py "https://github.com/niemasd/MultiVirusConsensus/raw/refs/heads/main/MultiVirusConsensus.py" && \
    chmod a+x /usr/local/bin/MultiVirusConsensus.py && \

    # clean up
    apt-get clean && rm -rf /root/.cache /tmp/*
