# OPUS dei 📊

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/magonji/opus-converter/issues)

**OPUS dei v2.0** — A powerful command-line tool for converting Bruker OPUS spectral files to accessible formats.

## 🎯 Overview

OPUS Converter is designed to efficiently process Bruker OPUS binary files and convert them into more accessible formats for spectral analysis. Whether you're dealing with dozens or tens of thousands of spectral files, this tool provides an intelligent, user-friendly interface with adaptive processing capabilities.

### Key Features

- 📁 **Recursive File Discovery** — Automatically finds OPUS files in all subdirectories
- ⚡ **Smart Batch Processing** — Adaptive interface that scales from small to massive datasets
- 🎯 **Flexible Output Formats** — Choose exactly the formats you need
- 📊 **Real-time Progress Tracking** — ETA calculations and processing statistics
- 🛡️ **Error Resilience** — Continues processing even when individual files fail
- 🎨 **Beautiful Terminal Interface** — Colourful, informative output using modern CLI design

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/magonji/opus-dei.git
cd opus-converter

# Install dependencies
pip install -r requirements.txt

# Run the converter
python opus_converter.py
```

### Requirements

- Python 3.7 or higher
- Required packages:
  ```
  numpy
  questionary
  colorama
  ```

## 📋 Usage

Simply run the script and follow the interactive prompts:

```bash
python opus_converter.py
```

The program will guide you through three main steps:

### 1. Directory Selection 📁
Browse and select the folder containing your OPUS files. The program recursively searches through all subdirectories to find files with numerical extensions (e.g., `spectrum.0`, `data.1`, `sample.15`).

### 2. Format Selection 🎯
Choose your preferred output format(s):

- **📊 Both .dpt and .mzz files** — Complete conversion with both formats
- **📈 Only .dpt files** — Full resolution tab-delimited text files
- **🗜️ Only .mzz files** — Space-efficient compressed format (1 cm⁻¹ resolution)

### 3. Batch Processing ⚡
The tool automatically adapts its interface based on your dataset size:
- **Small batches (≤100 files)**: Individual file progress tracking
- **Large batches (>100 files)**: Progress bar with ETA and statistics

## 📄 Output Formats

### .dpt Files (Data Point Table)
- **Format**: Tab-delimited text file
- **Content**: Full resolution wavenumber and absorption data
- **Use case**: Direct import into analysis software (Origin, MATLAB, Python, R)
- **File size**: Larger, but preserves all spectral information

### .mzz Files (Compressed Spectral Archive)
- **Format**: ZIP-compressed custom format
- **Content**: Rounded wavenumbers (1 cm⁻¹ resolution) with intensity data
- **Use case**: Long-term storage, large dataset management
- **File size**: Significantly smaller, optimised for storage efficiency

## 🧪 Technical Details

### OPUS File Structure
The program parses Bruker OPUS binary files by:
- Reading the file header to locate data blocks
- Extracting absorption spectra and metadata
- Processing wavenumber ranges and spectral parameters
- Handling various block types (sample, reference, parameters, etc.)

### Performance Optimisation
- **Memory efficient**: Processes files individually to handle large datasets
- **Error handling**: Robust error management with detailed reporting
- **Cross-platform**: Works on Windows, macOS, and Linux
- **Unicode support**: Handles international characters in file paths

## 🔧 Development

### Project Structure
```
opus-converter/
├── opus_converter.py          # Main application
├── requirements.txt           # Python dependencies
├── README.md                 # This file
└── examples/                 # Example files and usage
```

### Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📊 Use Cases

Perfect for:
- **Laboratory data archival** — Converting proprietary formats for long-term storage
- **Spectral database preparation** — Preparing datasets for machine learning or statistical analysis
- **Format migration** — Moving from OPUS to open formats for better accessibility
- **High-throughput processing** — Handling large batches from automated measurements
- **Cross-platform compatibility** — Accessing OPUS data on non-Windows systems

## ⚠️ Known Limitations

- Requires OPUS files to have numerical file extensions (standard OPUS convention)
- .mzz format uses 1 cm⁻¹ resolution (suitable for most FTIR applications)
- Designed specifically for absorption spectra (AB data blocks)

## 📞 Support

- **Issues**: Please report bugs or request features via [GitHub Issues](https://github.com/magonji/opus-converter/issues)
- **Email**: [mario.gonzalezjimenez@glasgow.ac.uk](mailto:mario.gonzalez-jimenez@glasgow.ac.uk)
- **Institution**: University of Glasgow

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🎓 Citation

If you use this software in your research, please consider citing:

```bibtex
@software{opus_converter_2025,
  author = {González-Jiménez, Mario},
  title = {OPUS Converter: A tool for converting Bruker OPUS spectral files},
  url = {https://github.com/magonji/opus-converter},
  version = {2.0},
  year = {2025},
  institution = {University of Glasgow}
}
```

## 🙏 Acknowledgments

- Developed at the University of Glasgow
- Inspired by the need for open-source spectral data processing tools
- Built with love for the spectroscopy community

---

**Made with ❤️ for the spectroscopy community**

*Transform your OPUS files with style* ✨
