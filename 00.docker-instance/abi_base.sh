!/usr/bin/env bash
# shellcheck shell=bash

# set -xv

## debugging, verbose mode
# Strict mode http://redsymbol.net/articles/unofficial-bash-strict-mode/
# set -euo pipefail
# IFS=$'\n\t'
# setopt NULL_GLOB

# Check if script is running on Amazon Linux or Ubuntu

LINUX_OS_TYPE=""
[[ $(uname -a) == *"amzn2023"* ]] &&
  echo "Running on Amazon Linux 2023" &&
  LINUX_OS_TYPE="amazon"

[[ $(uname -a) == *"Ubuntu"* ]] &&
  echo "Running on Ubuntu" &&
  LINUX_OS_TYPE="ubuntu" &&
  sudo ln -fs /usr/bin/apt-get /usr/local/bin/dnf

[[ $LINUX_OS_TYPE == "" ]] &&
  echo "Running on unknown OS" &&
  return 1

## ---- init - selinux - sshd

sudo dnf update -y
sudo dnf upgrade -y

sudo cat /etc/selinux/config |
  sed -E 's/^SELINUX=enforcing/SELINUX=permissive/g' >tmp
sudo mv tmp /etc/selinux/config

if [[ $(command -v setenforce) ]]; then
  sudo setenforce permissive
  sudo getenforce
fi

# https://stackoverflow.com/questions/8250379/sftp-on-linux-server-gives-error-received-message-too-long
sudo cp /etc/ssh/sshd_config /etc/ssh/sshd_config_orig &&
  sudo sed -i 's/#Port 22/Port 443/g' /etc/ssh/sshd_config &&
  sudo sed -i 's/^Subsystem.*/Subsystem sftp internal-sftp/g' /etc/ssh/sshd_config

sudo cat /etc/ssh/sshd_config | grep 443
sudo cat /etc/ssh/sshd_config | grep 'Subsystem'

sudo service sshd restart
sudo systemctl restart sshd

## ---- create required directories

mkdir -p "$HOME"/Dropbox &&
  mkdir -p "$HOME"/Dropbox/andrewrech &&
  mkdir -p "$HOME"/Dropbox/andrewrech/config &&
  mkdir -p "$HOME"/Dropbox/andrewrech/config/secrets &&
  mkdir -p "$HOME"/bin &&
  mkdir -p "$HOME"/.config &&
  mkdir -p "$HOME"/.config/rclone &&
  mkdir -p "$HOME"/.aws &&
  mkdir -p "$HOME"/Dropbox/andrewrech &&
  sudo mkdir -p /mnt &&
  sudo mkdir -p /mnt/scratch &&
  sudo mkdir -p /mnt/scratch-rclone-backup &&
  sudo mkdir -p /usr/share/fonts

sudo chown -R "$USER" /mnt/scratch

## --- required base packages

[[ $LINUX_OS_TYPE == "amazon" ]] &&
  sudo dnf groupinstall -y "Development Tools"

[[ $LINUX_OS_TYPE == "ubuntu" ]] &&
  sudo dnf install -y build-essential autoconf automake

sudo dnf install -y \
  awscli \
  git \
  parallel \
  python3-pip \
  tar \
  tmux \
  unzip \
  zstd \
  zsh &&
  sudo dnf clean all

[[ $LINUX_OS_TYPE == "amazon" ]] &&
  sudo dnf install -y \
    R-devel \
    R-java \
    python3-devel &&
  sudo dnf clean all

[[ $LINUX_OS_TYPE == "ubuntu" ]] &&
  sudo dnf install -y \
    r-base-dev \
    r-recommended \
    r-cran-r-java \
    python-all-dev &&
  sudo dnf clean all

pip3 install --upgrade pip
pip3 install \
  glances

## ---- cli helpers

[[ -d /tmp/autojump ]] && rm -rf /tmp/autojump
cd /tmp &&
  git clone 'https://github.com/wting/autojump' &&
  cd autojump &&
  ./install.py

[[ $LINUX_OS_TYPE == "amazon" ]] &&
  sudo dnf install -y \
    libX11-devel \
    xorg-x11-server-Xvfb

## ---- zsh

touch "$HOME"/.zshrc
sudo usermod --shell /usr/bin/zsh "$USER"
sudo usermod --shell /usr/bin/zsh root
sudo cp /usr/bin/zsh /usr/local/bin/zssudo yum install h

## ---- set locale

[[ $LINUX_OS_TYPE == "amazon" ]] &&
  sudo dnf install -y \
    glibc-locale-source \
    glibc-langpack-en &&
  sudo dnf clean all

[[ $LINUX_OS_TYPE == "ubuntu" ]] &&
  sudo dnf install -y \
    locales \
    language-pack-en

sudo localedef -v -c -i en_US -f UTF-8 en_US.UTF-8 || true

## ---- atuin

if ! command -v atuin >/dev/null 2>&1; then
  echo "Atuin not found in PATH. Installing Atuin..."
  curl --proto '=https' --tlsv1.2 -LsSf https://setup.atuin.sh | sh
fi

## ---- rclone latest requires to read config
if [[ $(uname -m) == "aarch64" ]]; then
  ARCH="arm64"
else
  ARCH="amd64"
fi

wget \
  "$(curl -s https://api.github.com/repos/rclone/rclone/releases/latest | grep -Po "https://github.com/rclone/rclone/releases/download/v[0-9\.]+/rclone-v[0-9\.]+-linux-$ARCH.zip")" &&
  unzip rclone-v*-linux-$ARCH.zip &&
  cd rclone-v*-linux-$ARCH &&
  sudo cp rclone /usr/local/bin &&
  cd ../ &&
  rm -rf rclone-v*-linux*

## ---- mosh

[[ $LINUX_OS_TYPE == "amazon" ]] &&
  sudo dnf install -y \
    protobuf-devel \
    ncurses-devel \
    openssl-devel

[[ $LINUX_OS_TYPE == "ubuntu" ]] &&
  sudo dnf install -y \
    libprotobuf-dev \
    libncurses-dev \
    libssl-dev

cd /tmp &&
  git clone "https://github.com/mobile-shell/mosh" &&
  cd mosh &&
  ./autogen.sh &&
  ./configure &&
  make -j$NPROC &&
  sudo make install &&
  cd /tmp &&
  sudo rm -rf /tmp/mosh

## ---- pzstd

cd /tmp &&
  git clone "https://github.com/facebook/zstd/" &&
  cd ./zstd/contrib/pzstd &&
  make -j$NRPOC &&
  sudo make install

## ---- docker

[[ $LINUX_OS_TYPE == "amazon" ]] &&
  sudo dnf install -y \
    docker

[[ $LINUX_OS_TYPE == "ubuntu" ]] &&
  sudo dnf install -y \
    docker.io \
    docker-buildx

sudo systemctl enable docker.service &&
  sudo service docker start &&
  sudo usermod -a -G docker "$USER" &&
  sudo usermod -a -G docker root &&
  sudo systemctl restart docker.service

## ---- NVIDIA driver, CUDA toolkit, and Container Toolkit

echo "Installing NVIDIA driver, CUDA toolkit, and container toolkit..."

dnf check-release-update
sudo dnf upgrade --releasever=latest -y
sudo dnf update -y
sudo dnf install -y dkms kernel-devel kernel-modules-extra
sudo systemctl enable dkms

cd /tmp
sudo dnf install -y vulkan-devel libglvnd-devel elfutils-libelf-devel xorg-x11-server-Xorg
# manually check for latest version
# https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=arm64-sbsa&Compilation=Native&Distribution=RHEL&target_version=9&target_type=runfile_local
wget https://developer.download.nvidia.com/compute/cuda/12.8.1/local_installers/cuda_12.8.1_570.124.06_linux_sbsa.run
chmod +x ./cuda*.run
sudo ./cuda_*.run --driver --toolkit --tmpdir=/var/tmp --silent

if (! dnf search nvidia | grep -q nvidia-container-toolkit); then
  sudo dnf config-manager --add-repo https://nvidia.github.io/libnvidia-container/stable/rpm/nvidia-container-toolkit.repo
fi
sudo dnf install --nogpgcheck -y nvidia-container-toolkit

sudo dnf install -y docker
sudo systemctl enable docker
sudo usermod -aG docker ec2-user

sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

sudo dnf clean all
sudo rm -rf /var/cache/yum
sudo reboot now
