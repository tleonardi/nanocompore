# Troubleshooting


### Missing libcufft.so

If you're getting error like `OSError: libcufft.so.11: cannot open shared object file: No such file or directory` it may be because your system doesn't have an NVidia GPU. A possible workaround is to use the CPU version of PyTorch:

```bash
$ pip3 install torch  --extra-index-url https://download.pytorch.org/whl/cpu
```

### Some transcripts fail with an Out of memory error

Please take a look at the [Common pitfalls section](/usage/#common-pitfalls).

### CUDA error: no kernel image is available for execution on the device

This error may happen when there's a mismatch between the GPU architecture and the CUDA toolkit that is installed. Please make sure that you have installed the correct driver version for your GPU model and the version of the CUDA toolkit is compatible with the drivers and the GPU model. You can also check stderr, since pytorch may output a more detailed error message there. For information on CUDA/GPU comptability refer to the [NVidia's official documentation](https://developer.nvidia.com/cuda/gpus).

