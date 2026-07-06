# OPUS dei 📊

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/magonji/opus-dei/issues)

**OPUS dei v3.0** — A powerful command-line tool for converting Bruker OPUS spectral files to accessible formats.

## 🎯 Overview

OPUS Converter is designed to efficiently process Bruker OPUS binary files and convert them into more accessible formats for spectral analysis. Whether you're dealing with dozens or tens of thousands of spectral files, this tool provides an intelligent, user-friendly interface with adaptive processing capabilities.

### Key Features

- 📁 **Recursive File Discovery** — Automatically finds OPUS files in all subdirectories
- ⚡ **Smart Batch Processing** — Adaptive interface that scales from small to massive datasets
- 🎯 **Flexible Output Formats** — Choose exactly the formats you need
- 🔎 **Metadata Inspector** — Read all the metadata stored in a single OPUS file on screen (grouped by block, with readable descriptions of common parameter codes) and optionally save it to a `.txt`
- 📊 **Real-time Progress Tracking** — ETA calculations and processing statistics
- 🛡️ **Error Resilience** — Continues processing even when individual files fail
- 🎨 **Full-screen TUI** — A smooth, retro terminal interface with a giant ASCII banner, arrow-key navigation, a filterable file browser and an in-place progress view (no scrolling)

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/magonji/opus-dei.git
cd opus-dei

# Install dependencies
pip install -r requirements.txt

# Run the converter
python opus-dei.py
```

### Requirements

- Python 3.7 or higher
- Required packages:
  ```
  numpy
  prompt_toolkit
  colorama
  ```

## 📋 Usage

Run the script and move through the full-screen interface with the keyboard
(`↑↓` to move, `⏎` to select, `←` to go back, `q` to quit). In the browser you
can navigate with the arrow keys, press `Tab` to type/paste a path directly,
and (when picking a file) press `/` to filter the list by text:

```bash
python opus-dei.py
```

On the first screen you choose what to do:

- **Convert spectra** — batch-convert a folder of OPUS files to `.dpt` / `.mzz`
  (the three steps below).
- **Inspect metadata** — pick a single OPUS file and read all of its metadata
  on screen, grouped by block with readable descriptions of the common
  parameter codes; press `s` to save it to a `.txt` next to the file.

### Convert spectra

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

``` bash
opus-dei/
├── opus-dei.py          # Main application
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

- **Issues**: Please report bugs or request features via [GitHub Issues](https://github.com/magonji/opus-dei/issues)
- **Email**: [mario.gonzalezjimenez@glasgow.ac.uk](mailto:mario.gonzalez-jimenez@glasgow.ac.uk)
- **Institution**: University of Glasgow

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🎓 Citation

If you use this software in your research, please consider citing:

```bibtex
@software{opus-dei_2025,
  author = {González-Jiménez, Mario},
  title = {OPUS dei: A tool for converting Bruker OPUS spectral files},
  url = {https://github.com/magonji/opus-dei},
  version = {3.0},
  year = {2025},
  institution = {University of Glasgow}
}
```

## 🙏 Acknowledgments

- Developed at the University of Glasgow
- Inspired by the need for open-source spectral data processing tools
- Built with love for the spectroscopy community

---

**Made with lots of ❤️ for the spectroscopy community/Bruker sufferers**

*Transform your OPUS files with style* ✨
