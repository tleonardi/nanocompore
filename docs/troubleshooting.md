# Troubleshooting


### Missing libcufft.so

If you're getting error like `OSError: libcufft.so.11: cannot open shared object file: No such file or directory` it may be because your system doesn't have an NVidia GPU. A possible workaround is to use the CPU version of PyTorch:

```bash
$ pip3 install torch  --extra-index-url https://download.pytorch.org/whl/cpu
```

