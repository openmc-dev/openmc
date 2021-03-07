#define catchCudaErrors(msg)                                                   \
  {                                                                            \
    __catchCudaErrors(msg, __FILE__, __LINE__);                                \
  }
inline void __catchCudaErrors(const char* message, const char* file, int line)
{
  auto error_code = cudaGetLastError();
  if (error_code != cudaSuccess) {
    fprintf(stderr, "%s: %s %s line %d\n", message,
      cudaGetErrorString(error_code), file, line);
    exit(error_code);
  }
}
