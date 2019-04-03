#include <stdlib.h>

double cfloat_endf(const char* buffer, int n)
{
  char arr[12]; // 11 characters plus a null terminator
  int j = 0; // current position in arr
  int found_decimal = 0;
  int found_exponent = 0;
  for (int i = 0; i < n; ++i) {
    // Skip whitespace characters
    char c = buffer[i];
    if (c == ' ') continue;

    if (found_decimal) {
      if (!found_exponent) {
        if (c == '+' || c == '-') {
          // In the case that we encounter +/- and we haven't yet encountered
          // e/E, we manually add it
          arr[j++] = 'e';
          found_exponent = 1;

        } else if (c == 'e' || c == 'E') {
          found_exponent = 1;
        }
      }
    } else if (c == '.') {
      found_decimal = 1;
    }

    // Copy character
    arr[j++] = c;
  }

  // Done copying. Add null terminator and convert to double
  arr[j] = '\0';
  return atof(arr);
}
