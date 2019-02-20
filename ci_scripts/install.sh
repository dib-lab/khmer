# Used to prepare OSX builds

if [ "${TRAVIS_OS_NAME}" != "osx" ]; then
  echo "This is not a OSX build"
  exit 1
fi

brew update;
brew install astyle;
brew install cppcheck;

pushd .
cd
mkdir -p download
cd download
wget https://repo.continuum.io/miniconda/Miniconda3-4.5.11-MacOSX-x86_64.sh \
  -O miniconda.sh
chmod +x miniconda.sh && ./miniconda.sh -b
cd ..
export PATH=/Users/travis/miniconda3/bin:$PATH
conda update --yes conda
popd

# Create a fresh environment
conda create -n testenv --yes python=3

source activate testenv

# Without this `make install-dependencies` fails
pip install --ignore-installed setuptools
