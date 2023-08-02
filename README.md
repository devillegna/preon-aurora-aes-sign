# preon-aurora-aes-sign

preon: aurora-aes based signature

## Prepare environments

1. Python packages

```
pip install -r requirements.txt
```

2. Pull submodules

```
git submodule update --init --recursive
```

3. Prepare shared library for FFT in C

```
make
```

4. Installing Jupyter Notebook(https://jupyter.org/) for running 'xxx.ipynb' testers. 

## Contents

- **aesR1CS** : R1CS for AES128.
- **aurora**  : The Aurora proving system.
- **fri_ldt** : FRI low degree test: a ZK proof of polynomial degrees.
- **commit**  : Merkel tree commit scheme.
- **libgf264fft**  : Python wrappers for C fft library.
- **polyeval** : Submodule pointing to C fft library.
- **preon**  : The Preon digital signature scheme.
- **utils**  : Some utilities: field arithmetic, hash functions, PRNGs, etc.


## Run examples

Please see/run preon-test.ipynb with jupyter notebook(https://jupyter.org/).

