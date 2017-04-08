# picojdec

Simple JPEG decoder implementation.

This project is intended for study and demonstration, not for practical use. The implementation is naive and inefficiency, so its decoding speed is VERY slow.


# Limitations
- support 'Baseline JPEG' only ('Progressive JPEG' is not supported)
- support Grayscale and YCbCr color format only (RGB/CYMK... and ICC profiles are not supported)
- support Huffman coding only (Arithmetic coding is not supported)
- output [PPM(portable pixmap)][ppm] image file only

[ppm]: https://en.wikipedia.org/wiki/Netpbm_format


# References
- [ITU-T(CCITT) Rec. T.81][t81] "Digital compression and coding of continuous-tone still images"
  - [free access version][t81-w3c] on www.w3.org site
- [ITU-T Rec. T.871][t871] "Digital compression and coding of continuous-tone still images: JPEG File Interchange Format (JFIF)"

[t81]: https://www.itu.int/rec/T-REC-T.81
[t81-w3c]: https://www.w3.org/Graphics/JPEG/itu-t81.pdf
[t871]: https://www.itu.int/rec/T-REC-T.871


# License
MIT License
