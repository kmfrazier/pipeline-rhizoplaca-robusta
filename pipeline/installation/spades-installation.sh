#IF NOT USING LINUX, FOLLOW INSTRUCTIONS HERE:
#https://ablab.github.io/spades/installation.html

if [[ "$OSTYPE" == "linux-gnu"* ]]; then 
    wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
    tar -xzf SPAdes-4.0.0-Linux.tar.gz
    mv SPAdes-4.0.0-Linux ../assembly/SPAdes-4.0.0-Linux
    cd ../assembly/SPAdes-4.0.0-Linux/bin/
    echo "SPAdes installation complete"
else 
    echo "Not linux; check file for installation instructions"
fi