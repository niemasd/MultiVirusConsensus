# Minimal Docker image for MultiVirusConsensus using Alpine base
FROM alpine:3.13.5
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# install MultiVirusConsensus and dependencies
RUN apk update && \
    apk add autoconf automake bash bzip2-dev curl-dev g++ gcc make musl-dev xz-dev zlib-dev && \

    # install htslib (needed for ViralConsensus)
    wget -qO- "https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2" | tar -xj && \
    cd htslib-* && autoreconf -i && ./configure && make && make install && cd .. && rm -rf htslib-* && \

    # install samtools
    wget -qO- "https://github.com/samtools/samtools/releases/download/1.23.1/samtools-1.23.1.tar.bz2" | tar -xj && \
    cd samtools-* && ./configure --without-curses && make && make install && cd .. && rm -rf samtools-* && \

    # install minimap2
    wget -qO- "https://github.com/lh3/minimap2/archive/refs/tags/v2.31.tar.gz" | tar -zx && \
    cd minimap2-* && make && chmod a+x minimap2 && mv minimap2 /usr/local/bin/minimap2 && cd .. && rm -rf minimap2-* && \

    # install ViralConsensus
    wget -qO- "https://github.com/niemasd/ViralConsensus/archive/refs/tags/1.0.1.tar.gz" | tar -zx && \
    cd ViralConsensus-* && make && mv viral_consensus /usr/local/bin/viral_consensus && cd .. && rm -rf ViralConsensus-* && \

    # install MultiVirusConsensus
    wget -q -O /usr/local/bin/MultiVirusConsensus.py "https://github.com/niemasd/MultiVirusConsensus/raw/refs/heads/main/MultiVirusConsensus.py" && \
    chmod a+x /usr/local/bin/MultiVirusConsensus.py && \

    # clean up
    apt-get clean && rm -rf /root/.cache /tmp/*
