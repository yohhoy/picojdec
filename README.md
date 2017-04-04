# picojdec

Simple JPEG decoder implementation.

This project is intended for study and demonstration, not for practical use. The implementation is naive and inefficiency, so its decoding speed is VERY slow.


# Limitations
- support 'Baseline JPEG' only ('Progressive JPEG' is not supported)
- support Grayscale and YCbCr color format only (RGB/CYMK... and ICC profiles not supported)
- support Huffman coding only (Arithmetic coding is not supported)


# References
- [ITU-T(CCITT) Rec. T.81][t81] "Digital compression and coding of continuous-tone still images"
- [ITU-T Rec. T.871][t871] "Digital compression and coding of continuous-tone still images: JPEG File Interchange Format (JFIF)"

[t81]: https://www.itu.int/rec/T-REC-T.81
[t871]: https://www.itu.int/rec/T-REC-T.871


# License
MIT License
