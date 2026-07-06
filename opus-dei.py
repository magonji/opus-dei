#!/usr/bin/env python3
"""
OPUS File Converter

This script converts Bruker OPUS spectral files to two formats:
1. .dpt files: Tab-delimited text files with wavenumber and absorption data
2. .mzz files: Compressed format with rounded wavenumbers for space efficiency

Author: Converted from a Jupyter notebook
Usage: python opus_converter.py
"""

import numpy as np
import struct
import os
import io
import zipfile
import threading
import contextlib
import datetime
from pathlib import Path
from colorama import init, Fore, Back, Style
import time

from prompt_toolkit.application import Application
from prompt_toolkit.key_binding import KeyBindings
from prompt_toolkit.layout import Layout
from prompt_toolkit.layout.containers import (
    HSplit, Window, ConditionalContainer, WindowAlign,
)
from prompt_toolkit.layout.controls import FormattedTextControl
from prompt_toolkit.layout.dimension import Dimension
from prompt_toolkit.styles import Style as PTStyle
from prompt_toolkit.widgets import TextArea
from prompt_toolkit.filters import has_focus, Condition

# Initialise colorama for cross-platform colour support
init(autoreset=True)


# --- Splash screen -------------------------------------------------------

# Giant "OPUS DEI" title (figlet 'ansi_shadow' font).
TITLE_ART = r"""
 ██████╗ ██████╗ ██╗   ██╗███████╗    ██████╗ ███████╗██╗
██╔═══██╗██╔══██╗██║   ██║██╔════╝    ██╔══██╗██╔════╝██║
██║   ██║██████╔╝██║   ██║███████╗    ██║  ██║█████╗  ██║
██║   ██║██╔═══╝ ██║   ██║╚════██║    ██║  ██║██╔══╝  ██║
╚██████╔╝██║     ╚██████╔╝███████║    ██████╔╝███████╗██║
 ╚═════╝ ╚═╝      ╚═════╝ ╚══════╝    ╚═════╝ ╚══════╝╚═╝
""".strip("\n").splitlines()


# Colour palette (pink -> gold gradient endpoints, plus accents).
GRAD_START = (255, 58, 107)
GRAD_END = (255, 170, 64)
ACCENT = "#ff3a6b"
GOLD = "#ffaa40"
DIM = "#8a8a96"


def _hex(rgb):
    """Convert an (r, g, b) tuple to a #rrggbb string."""
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def title_fragments(pad="   "):
    """Return the giant OPUS DEI title as prompt_toolkit (style, text) fragments.

    Each line gets its own colour along a vertical pink -> gold gradient.
    prompt_toolkit downgrades the colours automatically on terminals with a
    limited palette, so no manual fallback is needed.
    """
    fragments = []
    n = max(len(TITLE_ART) - 1, 1)
    for i, line in enumerate(TITLE_ART):
        r = round(GRAD_START[0] + (GRAD_END[0] - GRAD_START[0]) * i / n)
        g = round(GRAD_START[1] + (GRAD_END[1] - GRAD_START[1]) * i / n)
        b = round(GRAD_START[2] + (GRAD_END[2] - GRAD_START[2]) * i / n)
        fragments.append((f"fg:{_hex((r, g, b))} bold", pad + line + "\n"))
    return fragments


class OpusFileReader(dict):
    """
    A class to read and parse Bruker OPUS spectral files.
    
    This class inherits from dict to store the parsed data blocks
    and provides methods to read the binary OPUS file format.
    """
    
    def __init__(self, filepath):
        """
        Initialise the OPUS file reader.
        
        Args:
            filepath (str): Path to the OPUS file to read
        """
        super().__init__()
        
        with open(filepath, 'rb') as opus_file:
            self.raw_data = opus_file.read()
        
        self.total_data_length = len(self.raw_data)
        self._read_file_header()
        
        self.data_blocks = []
        self.parameter_list = []

    def _read_file_header(self):
        """
        Read and parse the OPUS file header to extract data block information.
        
        The header contains information about data blocks including their
        offsets, sizes, types, and channels.
        """
        header_size = 504
        self.header = self.raw_data[0:header_size]

        # Initialise lists to store header information
        self.block_offsets = []
        self.chunk_sizes = []
        self.block_types = []
        self.channel_types = []
        self.text_types = []

        cursor_position = 32
        
        while cursor_position > 0:
            start_index = cursor_position
            end_index = start_index + 4

            if end_index <= header_size:
                # Read block offset (4 bytes, little-endian unsigned int)
                block_offset = struct.unpack('<I', self.header[start_index:end_index])[0]
                
                if block_offset > 0:
                    self.block_offsets.append(block_offset)
                    
                    # Read chunk size (4 bytes before offset)
                    size_start = cursor_position - 4
                    size_end = size_start + 4
                    chunk_size = struct.unpack('<I', self.header[size_start:size_end])[0]
                    self.chunk_sizes.append(chunk_size)
                    
                    # Read data type (1 byte, 8 bytes before offset)
                    type_start = cursor_position - 8
                    type_end = type_start + 1
                    data_type = struct.unpack('<B', self.header[type_start:type_end])[0]
                    self.block_types.append(data_type)

                    # Read channel type (1 byte, 7 bytes before offset)
                    channel_start = cursor_position - 7
                    channel_end = channel_start + 1
                    channel_type = struct.unpack('<B', self.header[channel_start:channel_end])[0]
                    self.channel_types.append(channel_type)

                    # Read text type (1 byte, 6 bytes before offset)
                    text_start = cursor_position - 6
                    text_end = text_start + 1
                    text_type = struct.unpack('<B', self.header[text_start:text_end])[0]
                    self.text_types.append(text_type)

                    next_offset = block_offset + 4 * chunk_size
                    
                    if next_offset >= self.total_data_length:
                        # Next offset would reach end of file
                        cursor_position = -1
                    else:
                        cursor_position += 12
                else:
                    cursor_position = -1
            else:
                cursor_position = -1

    def read_all_data_blocks(self):
        """
        Read and process all data blocks found in the OPUS file.
        
        This method processes each data block according to its type and
        stores the results in the dictionary with descriptive names.
        """
        num_blocks = len(self.block_offsets)
        
        for block_index in range(num_blocks):
            raw_chunk = self._read_raw_chunk(block_index)
            chunk_size = self.chunk_sizes[block_index]
            block_type = self.block_types[block_index]
            text_type = self.text_types[block_index]
            channel_type = self.channel_types[block_index]
            
            data_block = DataBlock(
                raw_chunk=raw_chunk, 
                chunk_size=chunk_size,
                block_type=block_type, 
                text_type=text_type
            )
            
            self.data_blocks.append(data_block)
            block_name = self._determine_block_name(block_type, text_type, channel_type)
            
            if block_name:
                self[block_name] = data_block
                
                # Create parameter entry for the block
                parameter = {
                    'name': block_name, 
                    'type': 'group',
                    'children': data_block.parameter_list
                }
                self.parameter_list.append(parameter)

        # Generate wavenumber axis if absorption data is available
        if 'AB Data Parameter' in self.keys():
            first_wavenumber = self['AB Data Parameter']['FXV']
            last_wavenumber = self['AB Data Parameter']['LXV']
            num_points = self['AB Data Parameter']['NPT']
            self['WN'] = np.linspace(first_wavenumber, last_wavenumber, num_points)

    def _determine_block_name(self, block_type, text_type, channel_type):
        """
        Determine the descriptive name for a data block based on its type codes.
        
        Args:
            block_type (int): The block type identifier
            text_type (int): The text type identifier  
            channel_type (int): The channel type identifier
            
        Returns:
            str: Descriptive name for the block, or None if not recognised
        """
        if block_type == 0:
            # Text information blocks
            text_block_names = {
                8: 'Info Block',
                104: 'History',
                152: 'Curve Fit',
                168: 'Signature',
                240: 'Integration Method'
            }
            return text_block_names.get(text_type, 'Text Information')
            
        elif block_type == 7:
            # Single channel sample spectra
            channel_names = {4: 'ScSm', 8: 'IgSm', 12: 'PhSm'}
            if channel_type in channel_names:
                self[channel_names[channel_type]] = np.array(self.data_blocks[-1].spectral_values)
            return None
            
        elif block_type == 11:
            # Single channel reference spectra
            channel_names = {4: 'ScRf', 8: 'IgRf'}
            if channel_type in channel_names:
                self[channel_names[channel_type]] = np.array(self.data_blocks[-1].spectral_values)
            return None
            
        elif block_type == 15:
            # Absorption spectrum
            self['AB'] = np.array(self.data_blocks[-1].spectral_values)
            return None
            
        elif block_type == 23:
            # Sample data parameters
            parameter_names = {
                4: 'ScSm Data Parameter',
                8: 'IgSm Data Parameter',
                12: 'PhSm Data Parameter'
            }
            return parameter_names.get(channel_type)
            
        elif block_type == 27:
            # Reference data parameters
            parameter_names = {4: 'ScRf Data Parameter', 8: 'IgRf Data Parameter'}
            return parameter_names.get(channel_type)
            
        elif block_type == 31:
            return 'AB Data Parameter'
        elif block_type == 32:
            return 'Instrument'
        elif block_type == 40:
            return 'Instrument (Rf)'
        elif block_type == 48:
            return 'Acquisition'
        elif block_type == 56:
            return 'Acquisition (Rf)'
        elif block_type == 64:
            return 'Fourier Transformation'
        elif block_type == 72:
            return 'Fourier Transformation (Rf)'
        elif block_type == 96:
            return 'Optik'
        elif block_type == 104:
            return 'Optik (Rf)'
        elif block_type == 160:
            return 'Sample'
        else:
            print(f"{Fore.YELLOW}Warning: Unknown block type {block_type}{Style.RESET_ALL}")
            return None
    
    def _read_raw_chunk(self, block_index):
        """
        Extract raw data chunk for a specific block.
        
        Args:
            block_index (int): Index of the block to read
            
        Returns:
            bytes: Raw data chunk
        """
        start_pos = self.block_offsets[block_index]
        end_pos = start_pos + 4 * self.chunk_sizes[block_index]
        return self.raw_data[start_pos:end_pos]


class DataBlock(dict):
    """
    Represents a single data block from an OPUS file.
    
    This class handles parsing of different types of data blocks
    including spectral data, parameters, and text information.
    """
    
    def __init__(self, **kwargs):
        """
        Initialise a data block.
        
        Args:
            raw_chunk (bytes): Raw binary data for this block
            chunk_size (int): Size of the chunk in 4-byte units
            block_type (int): Type identifier for this block
            text_type (int): Text type identifier (default: -1)
        """
        super().__init__()
        
        self.text_type = kwargs.get('text_type', -1)
        self.raw_chunk = kwargs.get('raw_chunk')
        self.chunk_size = kwargs.get('chunk_size')
        self.block_type = kwargs.get('block_type')
        
        self.parameter_list = []
        self.spectral_values = None
        self.text_content = None
        
        self._parse_chunk_data()

    def _parse_chunk_data(self):
        """
        Parse the raw chunk data based on the block type.
        """
        if self.block_type == 0:
            if self.text_type == 8:
                # Info block with parameters
                self._parse_parameters()
            else:
                # Text content (history, etc.)
                self._parse_text_content()
        elif self.block_type in [7, 11, 15]:
            # Spectral data blocks
            self._parse_spectral_data()
        elif self.block_type in [23, 27, 31, 32, 40, 48, 64, 96, 104, 160]:
            # Parameter blocks
            self._parse_parameters()
        else:
            # Default to parameter parsing
            self._parse_parameters()
    
    def _parse_parameters(self):
        """
        Parse parameter data from the chunk.
        
        Parameters are stored as name-value pairs with type information.
        """
        cursor = 0
        parameter_types = ['int', 'float', 'str', 'str', 'str']

        while cursor >= 0:
            # Read parameter name (3 bytes)
            name_start = cursor
            name_end = name_start + 3

            try:
                param_name = self.raw_chunk[name_start:name_end].decode("utf-8")
            except UnicodeDecodeError:
                print(f"{Fore.YELLOW}Warning: Error decoding parameter name{Style.RESET_ALL}")
                break

            if param_name == 'END':
                break

            # Read parameter type (2 bytes, little-endian)
            type_start = cursor + 4
            type_end = type_start + 2
            type_index = struct.unpack('<H', self.raw_chunk[type_start:type_end])[0]

            try:
                param_type = parameter_types[type_index]
            except IndexError:
                print(f"{Fore.YELLOW}Warning: Unknown parameter type index {type_index}{Style.RESET_ALL}")
                break

            # Read parameter size (2 bytes, little-endian)
            size_start = cursor + 6
            size_end = size_start + 2
            param_size = struct.unpack('<H', self.raw_chunk[size_start:size_end])[0]

            # Read parameter value
            value_start = cursor + 8
            value_end = value_start + 2 * param_size
            raw_value = self.raw_chunk[value_start:value_end]

            # Parse value based on type
            param_value = self._parse_parameter_value(type_index, raw_value)
            
            # Store parameter
            self[param_name] = param_value
            
            parameter_info = {
                'name': param_name,
                'value': param_value,
                'type': param_type
            }
            self.parameter_list.append(parameter_info)

            cursor = cursor + 8 + 2 * param_size

    def _parse_parameter_value(self, type_index, raw_value):
        """
        Parse a parameter value based on its type index.
        
        Args:
            type_index (int): Type identifier for the parameter
            raw_value (bytes): Raw binary value data
            
        Returns:
            Parsed parameter value (int, float, or str)
        """
        if type_index == 0:
            # Integer (4 bytes, little-endian)
            return struct.unpack('<i', raw_value)[0]
        elif type_index == 1:
            # Double precision float (8 bytes, little-endian)
            return struct.unpack('<d', raw_value)[0]
        elif type_index in [2, 3, 4]:
            # String types (null-terminated, latin-1 encoding)
            null_pos = raw_value.find(b'\x00')
            if null_pos != -1:
                return raw_value[:null_pos].decode("latin-1")
            else:
                return raw_value.decode("latin-1")
        else:
            # Unknown type, return raw bytes
            return raw_value

    def _parse_spectral_data(self):
        """
        Parse spectral data as an array of floats.
        
        Spectral data is stored as little-endian 32-bit floats.
        """
        format_string = '<' + str(self.chunk_size) + 'f'
        self.spectral_values = struct.unpack(format_string, self.raw_chunk)

    def _parse_text_content(self):
        """
        Parse text content using latin-1 encoding.
        """
        self.text_content = self.raw_chunk.decode('latin-1')


def find_opus_files(directory_path):
    """
    Find all OPUS files in the specified directory and subdirectories.
    
    OPUS files are identified by having a filename that ends with a digit.
    
    Args:
        directory_path (str): Path to the directory to search
        
    Returns:
        list: List of paths to OPUS files found
    """
    opus_files = []
    search_path = Path(directory_path)
    
    for file_path in search_path.rglob('*'):
        if file_path.is_file() and file_path.name[-1].isdigit():
            opus_files.append(str(file_path))
    
    return opus_files


def convert_opus_file(opus_filepath, output_formats, show_individual_files=False):
    """
    Convert a single OPUS file to the specified formats.
    
    Args:
        opus_filepath (str): Path to the OPUS file to convert
        output_formats (list): List of formats to generate ('dpt', 'mzz', or both)
        show_individual_files (bool): Whether to show individual file processing messages
        
    Returns:
        bool: True if conversion was successful, False otherwise
    """
    try:
        if show_individual_files:
            print(f"{Fore.CYAN}📄 Processing:{Style.RESET_ALL} {Path(opus_filepath).name}")
        
        # Read and parse the OPUS file
        opus_reader = OpusFileReader(opus_filepath)
        opus_reader.read_all_data_blocks()

        # Extract absorption spectrum data
        absorption_spectrum = opus_reader["AB"]

        # Extract wavenumber range parameters
        first_wavenumber = opus_reader["AB Data Parameter"]["FXV"]
        last_wavenumber = opus_reader["AB Data Parameter"]["LXV"]

        # Build the wavenumber axis with linspace: exactly one point per data
        # value and exact endpoints. (The previous np.arange with a float step
        # drifted at the last point and, for ~1 in 4 files, produced one extra
        # point — which then made column_stack fail, silently dropping those
        # files from the conversion.)
        wavenumbers = np.linspace(
            first_wavenumber, last_wavenumber, len(absorption_spectrum)
        )

        # Create the full spectrum array (wavenumber, absorption)
        full_spectrum = np.column_stack((wavenumbers, absorption_spectrum))

        # Save the full spectrum to a .dpt file if requested
        if 'dpt' in output_formats:
            dpt_filepath = opus_filepath + ".dpt"
            np.savetxt(dpt_filepath, full_spectrum, fmt="%10.5f", delimiter="\t")
            if show_individual_files:
                print(f"  {Fore.GREEN}✓{Style.RESET_ALL} Created .dpt file")

        # Create compressed .mzz format if requested
        if 'mzz' in output_formats:
            # Create compressed version with rounded wavenumbers (1 cm⁻¹ resolution)
            rounded_spectrum = [[int(round(point[0])), point[1]] for point in full_spectrum]

            # Remove duplicate wavenumbers, keeping only the last occurrence
            compressed_spectrum = []
            for i in range(len(rounded_spectrum) - 1):
                if rounded_spectrum[i][0] != rounded_spectrum[i + 1][0]:
                    compressed_spectrum.append(rounded_spectrum[i])
            compressed_spectrum.append(rounded_spectrum[-1])  # Add the final point

            # Prepare data for export: start wavenumber, end wavenumber, number of points, then intensities
            export_data = [
                compressed_spectrum[0][0],      # First wavenumber
                compressed_spectrum[-1][0],     # Last wavenumber  
                len(compressed_spectrum)        # Number of data points
            ]
            export_data.extend([round(point[1], 4) for point in compressed_spectrum])

            # Save compressed data to temporary file
            temp_filepath = opus_filepath + ".tmp"
            np.savetxt(temp_filepath, export_data, fmt="%10.4f", delimiter="\t")

            # Create zip-compressed .mzz file
            mzz_filepath = opus_filepath + ".mzz"
            with zipfile.ZipFile(mzz_filepath, mode='w', compression=zipfile.ZIP_DEFLATED) as zip_file:
                zip_file.write(temp_filepath, arcname=Path(temp_filepath).name)

            # Clean up temporary file
            os.remove(temp_filepath)
            if show_individual_files:
                print(f"  {Fore.GREEN}✓{Style.RESET_ALL} Created .mzz file")
        
        return True

    except Exception as error:
        if show_individual_files:
            print(f"  {Fore.RED}✗ Error:{Style.RESET_ALL} {error}")
        return False


def _fmt_duration(seconds):
    """Format a duration in seconds as a compact human string."""
    if seconds is None:
        return "—"
    if seconds >= 60:
        return f"{int(seconds // 60)}m {int(seconds % 60)}s"
    if seconds >= 1:
        return f"{seconds:.1f}s"
    return f"{seconds * 1000:.0f}ms"


# What the user can do, chosen right after the welcome splash.
MODE_OPTIONS = [
    ("Convert spectra", "export .dpt / .mzz files from a folder", "convert"),
    ("Inspect metadata", "read all metadata from a single OPUS file", "inspect"),
]

# Output format options presented in the format-selection step.
FORMAT_OPTIONS = [
    ("Both .dpt and .mzz", "full resolution + compressed archive", ["dpt", "mzz"]),
    ("Only .dpt", "full resolution tab-delimited text", ["dpt"]),
    ("Only .mzz", "compressed 1 cm⁻¹ ZIP archive", ["mzz"]),
]

SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"

# Human-readable descriptions for the most common Bruker OPUS parameter codes.
# Codes not listed here are shown as-is (the raw three-letter code).
PARAM_DESCRIPTIONS = {
    # Data / axis parameters
    "FXV": "First X value (wavenumber)",
    "LXV": "Last X value (wavenumber)",
    "NPT": "Number of data points",
    "MXY": "Maximum Y value",
    "MNY": "Minimum Y value",
    "DXU": "X units",
    "DPF": "Data point format",
    "CSF": "Y-scaling factor",
    "DAT": "Measurement date",
    "TIM": "Measurement time",
    "DER": "Derivative",
    # Acquisition
    "RES": "Resolution (cm⁻¹)",
    "NSS": "Number of sample scans",
    "NSR": "Number of background scans",
    "AQM": "Acquisition mode",
    "PHR": "Phase resolution",
    "PHZ": "Phase correction mode",
    "APF": "Apodisation function",
    "ZFF": "Zero-filling factor",
    "DLY": "Stabilisation delay (s)",
    "HFW": "Wanted high-frequency limit",
    "LFW": "Wanted low-frequency limit",
    "GFW": "Number of good forward scans",
    "GBW": "Number of good backward scans",
    "BFW": "Number of bad forward scans",
    "BBW": "Number of bad backward scans",
    "PGN": "Preamplifier gain",
    "SGN": "Signal gain (sample)",
    "RGN": "Signal gain (reference)",
    # Fourier transform
    "FTS": "Number of FT points",
    "NLI": "Nonlinearity correction",
    # Optics
    "APT": "Aperture setting",
    "BMS": "Beamsplitter",
    "DTC": "Detector",
    "SRC": "Source",
    "VEL": "Scanner velocity",
    "OPF": "Optical filter",
    "PGR": "Preamplifier gain (optics)",
    "HPF": "High-pass filter",
    "LPF": "Low-pass filter",
    "CHN": "Measurement channel",
    # Instrument
    "INS": "Instrument type",
    "HFL": "High folding limit",
    "LFL": "Low folding limit",
    "LWN": "Laser wavenumber",
    "ASG": "Absolute signal gain",
    "ARG": "Absolute reference gain",
    "PKA": "Peak amplitude",
    "PKL": "Peak location",
    "SSP": "Sample spacing divisor",
    "RSN": "Running sample number",
    "GRB": "Reference/background number",
    # Sample / operator
    "SNM": "Sample name",
    "SFM": "Sample form",
    "CNM": "Operator name",
    "HIS": "History",
    "XPP": "Experiment path",
    "EXP": "Experiment name",
    "BLD": "Building",
    "CPY": "Company",
    "DPM": "Department",
    "LCT": "Location",
    # Accessory / channels
    "ACC": "Accessory",
    "RCH": "Reference channel",
    # Result / frequency range
    "PLF": "Result spectrum type",
    "HFQ": "Frequency limit, high (cm⁻¹)",
    "LFQ": "Frequency limit, low (cm⁻¹)",
    "DUR": "Measurement duration (s)",
    "DEL": "Delay before measurement (s)",
    # Instrument status / environment
    "SRN": "Instrument serial number",
    "VSN": "Firmware version",
    "IST": "Instrument status",
    "HUM": "Relative humidity (%)",
    "TSC": "Scanner temperature (°C)",
    "SRT": "Measurement start time",
}


def describe_param(code):
    """Return a human-readable description for an OPUS parameter code, or None."""
    return PARAM_DESCRIPTIONS.get(code.strip().upper())


# Parameter codes whose value is a Unix timestamp (seconds since the epoch).
_TIMESTAMP_CODES = {"SRT"}


def format_param_value(code, value):
    """Present a parameter value for display (e.g. timestamps as dates)."""
    if code.strip().upper() in _TIMESTAMP_CODES:
        try:
            dt = datetime.datetime.fromtimestamp(float(value))
            return dt.strftime("%Y-%m-%d %H:%M:%S")
        except (ValueError, OverflowError, OSError):
            pass
    return value


def _looks_like_text(value):
    """True if a string is mostly printable (so it isn't binary noise)."""
    if not value:
        return False
    printable = sum(1 for c in value if c.isprintable() or c in "\n\r\t")
    return printable / len(value) > 0.85


def _clean_value(value):
    """Strip NUL padding and control characters from a parameter value."""
    if not isinstance(value, str):
        return value
    return "".join(c for c in value if c.isprintable() or c in "\t\n").strip()


def metadata_text_lines(value):
    """Split a text-block value into readable lines (tabs become line breaks)."""
    return [line for line in str(value).replace("\t", "\n").splitlines() if line.strip()]


def collect_metadata(opus_filepath):
    """Read an OPUS file and return its metadata grouped by block.

    Returns a list of (block_name, [(code, value), …]) tuples, ordered as they
    appear in the file. Raises on unreadable files.
    """
    reader = OpusFileReader(opus_filepath)
    reader.read_all_data_blocks()

    blocks = []
    for group in reader.parameter_list:
        name = group["name"]
        params = [(p["name"], _clean_value(p["value"])) for p in group["children"]]
        if not params:
            # Text-only block (e.g. History): expose its text content, but only
            # if it really is text — some blocks decode to binary noise.
            data_block = reader.get(name)
            text = getattr(data_block, "text_content", None)
            if text and text.strip() and _looks_like_text(text):
                cleaned = _clean_value(text)
                if cleaned:
                    params = [("", cleaned)]
        if params:
            blocks.append((name, params))
    return blocks


def format_metadata_text(opus_filepath, blocks):
    """Render collected metadata as a human-readable plain-text report."""
    lines = [
        f"OPUS metadata — {Path(opus_filepath).name}",
        f"Source: {opus_filepath}",
        "=" * 64,
    ]
    for name, params in blocks:
        lines.append("")
        lines.append(f"[{name}]")
        for code, value in params:
            if not code:  # text block
                for text_line in metadata_text_lines(value):
                    lines.append(f"    {text_line}")
                continue
            desc = describe_param(code)
            label = f"{code}  ({desc})" if desc else code
            lines.append(f"  {label} = {format_param_value(code, value)}")
    return "\n".join(lines) + "\n"


class OpusDeiApp:
    """A full-screen, navigable TUI for the OPUS dei converter.

    The interface is a fixed screen (banner on top, a body that changes with
    the current step, and a status bar with shortcuts at the bottom). Nothing
    scrolls — the whole screen repaints in place, giving a smooth, modern feel.
    The spectral parsing/conversion engine is untouched and runs in a
    background thread so the UI stays responsive.
    """

    def __init__(self):
        # welcome | mode | folder | format | confirm | progress | summary | metadata
        self.step = "welcome"
        self.mode = "convert"          # convert | inspect
        self.mode_index = 0

        # Folder / file browser state.
        self.current_dir = Path.cwd()
        self.entries = []
        self.entry_index = 0

        # Format selection state.
        self.format_index = 0

        # Selections.
        self.selected_folder = None
        self.selected_formats = None

        # Direct path-entry state (folder step).
        self.path_error = None

        # Metadata inspector state.
        self.selected_file = None
        self.metadata = None           # list of (block, [(code, value), …])
        self.meta_error = None
        self.meta_scroll = 0
        self.meta_saved = None         # path of the saved .txt (confirmation)

        # Scanning / conversion state (shared with worker threads).
        self.opus_files = []
        self.scanning = False
        self.converting = False
        self.done_count = 0
        self.success = 0
        self.failed = 0
        self.start_time = None
        self.total_time = None

        self.app = self._build_application()

    # -- Terminal / helpers ------------------------------------------------

    @property
    def _spinner(self):
        return SPINNER[int(time.time() * 10) % len(SPINNER)]

    def _build_entries(self):
        """List selectable items for the browser at the current path.

        In convert mode this is a folder picker; in inspect mode it also lists
        the OPUS files (names ending in a digit) so one can be selected.
        """
        if self.mode == "convert":
            entries = [
                ("select", "Use this folder", self.current_dir),
                ("type", "Type or paste a folder path…", None),
            ]
        else:
            entries = [("type", "Type or paste a file path…", None)]

        parent = self.current_dir.parent
        if parent != self.current_dir:
            entries.append(("parent", ".. (parent folder)", parent))
        try:
            children = list(self.current_dir.iterdir())
        except (PermissionError, OSError):
            children = []

        for p in sorted((c for c in children if c.is_dir()), key=lambda p: p.name.lower()):
            entries.append(("dir", p.name + "/", p))

        if self.mode == "inspect":
            opus = sorted(
                (c for c in children if c.is_file() and c.name[-1:].isdigit()),
                key=lambda p: p.name.lower(),
            )
            for p in opus:
                entries.append(("file", p.name, p))

        self.entries = entries
        self.entry_index = 0

    def _visible_entries(self):
        """Entries after applying the live text filter (items only, not actions)."""
        text = self.filter_field.text.strip().lower()
        if not text:
            return self.entries
        visible = []
        for entry in self.entries:
            kind = entry[0]
            if kind in ("select", "type", "parent"):
                visible.append(entry)
            elif text in entry[1].lower():
                visible.append(entry)
        return visible

    # -- Rendering ---------------------------------------------------------

    def render_title(self):
        # No left padding: the window centres the banner horizontally.
        frags = title_fragments(pad="")
        # Drop the final newline so the art fits an exact-height window.
        style, text = frags[-1]
        frags[-1] = (style, text.rstrip("\n"))
        return frags

    def render_subtitle(self):
        return [
            (f"fg:{GOLD}", "Bruker OPUS spectral converter"),
            (f"fg:{DIM}", "   ·   "),
            (f"fg:{GOLD}", "v4.0"),
            (f"fg:{DIM}", "   ·   Mario González-Jiménez · University of Glasgow"),
        ]

    def render_body(self):
        return getattr(self, f"_body_{self.step}")()

    def _body_welcome(self):
        # The welcome splash is drawn by dedicated centred windows (banner,
        # subtitle, hint) so the whole thing can be centred; the body is empty.
        return []

    def render_welcome_hint(self):
        return [(f"fg:{GOLD} bold", "Press ⏎ to begin")]

    def _welcome_top_height(self):
        """Exact top padding that centres the welcome splash (0 elsewhere).

        The greedy bottom filler takes the matching space below, so the splash
        block lands dead centre and the status bar stays on the bottom edge.
        """
        if self.step != "welcome":
            return Dimension.exact(0)
        try:
            rows = self.app.output.get_size().rows
        except Exception:
            return Dimension.exact(0)
        block = len(TITLE_ART) + 1 + 1 + 1  # banner + subtitle + gap + hint
        return Dimension.exact(max((rows - block - 1) // 2, 0))  # -1 status bar

    def _is_typing_path(self):
        """True while the direct path-entry field holds keyboard focus."""
        try:
            return self.app.layout.has_focus(self.path_field)
        except Exception:
            return False

    def _is_filtering(self):
        """True while the list-filter field holds keyboard focus."""
        try:
            return self.app.layout.has_focus(self.filter_field)
        except Exception:
            return False

    def _focus_path_input(self):
        """Move focus to the path field for direct entry (empty, ready to type)."""
        self.path_error = None
        self.path_field.text = ""  # start empty so a pasted path just works
        self.app.layout.focus(self.path_field)

    def _accept_path(self, buffer):
        """Handle a path typed into the direct-entry field (Enter pressed)."""
        raw = buffer.text.strip()
        if not raw:
            return False
        candidate = Path(raw).expanduser()
        try:
            resolved = candidate.resolve()
        except (OSError, RuntimeError):
            resolved = candidate

        if resolved.is_file():
            self.path_error = None
            self.app.layout.focus(self.body_window)
            if self.mode == "inspect":
                self._open_file(resolved)          # open the file directly
            else:
                self.current_dir = resolved.parent  # convert: use its folder
                self._build_entries()
            return False
        if resolved.is_dir():
            self.current_dir = resolved
            self._build_entries()
            self.path_error = None
            self.app.layout.focus(self.body_window)
            return False

        self.path_error = f"Not found: {raw}"
        return True  # keep the text so the user can fix it

    def _body_mode(self):
        out = [("", "\n")]
        out.append((f"fg:{GOLD} bold", "   What would you like to do?\n\n"))
        for i, (name, desc, _value) in enumerate(MODE_OPTIONS):
            selected = i == self.mode_index
            if selected:
                out.append((f"fg:#ffffff bg:{ACCENT} bold", f"   › {name}"))
                out.append((f"fg:#ffffff bg:{ACCENT}", " " * max(20 - len(name), 1)))
                out.append((f"fg:#ffffff bg:{ACCENT}", f"{desc}"))
                out.append((f"fg:#ffffff bg:{ACCENT}", " " * 2 + "\n"))
            else:
                out.append(("", f"     {name}"))
                out.append(("", " " * max(20 - len(name), 1)))
                out.append((f"fg:{DIM}", f"{desc}\n"))
        return out

    def _body_folder(self):
        inspect = self.mode == "inspect"
        target = "OPUS file" if inspect else "folder"
        out = [("", "\n")]
        heading = ("Pick the OPUS file to inspect"
                   if inspect else "Choose the folder that contains your OPUS files")
        out.append((f"fg:{GOLD} bold", f"   {heading}\n"))

        if self._is_typing_path():
            out.append((f"fg:{DIM}", "   Type or paste the full path below, then press "))
            out.append((f"fg:{GOLD}", "⏎"))
            out.append((f"fg:{DIM}", "  (Esc to cancel)\n"))
            example = str(Path.home() / "Documents" / "OPUS data")
            out.append((f"fg:{DIM}", f"   e.g.  {example}\n\n"))
        else:
            out.append((f"fg:{DIM}", "   Browse with "))
            out.append((f"fg:{GOLD}", "↑↓"))
            out.append((f"fg:{DIM}", " and "))
            out.append((f"fg:{GOLD}", "⏎"))
            if inspect:
                out.append((f"fg:{DIM}", ", press "))
                out.append((f"fg:{GOLD}", "/"))
                out.append((f"fg:{DIM}", " to filter, "))
                out.append((f"fg:{GOLD}", "Tab"))
                out.append((f"fg:{DIM}", " to type a path.\n\n"))
            else:
                out.append((f"fg:{DIM}", ", or choose "))
                out.append((f"fg:{GOLD}", "“Type or paste a folder path”"))
                out.append((f"fg:{DIM}", ".\n\n"))

        out.append((f"fg:{DIM}", "   Currently in:  "))
        out.append(("", str(self.current_dir) + "\n"))
        if self.path_error:
            out.append((f"fg:{ACCENT}", "   ✗ " + self.path_error + "\n"))
        out.append(("", "\n"))

        entries = self._visible_entries()
        if inspect and not any(e[0] == "file" for e in entries):
            note = ("no OPUS files match the filter"
                    if self.filter_field.text.strip() else "no OPUS files in this folder")
            out.append((f"fg:{DIM}", f"   ({note})\n"))

        # Simple viewport so long listings never overflow the screen.
        max_visible = 12
        total = len(entries)
        if self.entry_index >= total:
            self.entry_index = max(total - 1, 0)
        start = max(0, min(self.entry_index - max_visible // 2, total - max_visible))
        start = max(start, 0)
        visible = entries[start:start + max_visible]

        if start > 0:
            out.append((f"fg:{DIM}", "     ↑ …\n"))
        icons = {"select": "✓ ", "type": "⌨  ", "parent": "↩ ", "dir": "📁 ", "file": "📄 "}
        for i, (kind, label, _path) in enumerate(visible, start=start):
            selected = i == self.entry_index
            text = f"{icons[kind]}{label}"
            if selected:
                out.append((f"fg:#ffffff bg:{ACCENT} bold", f"   › {text}"))
                out.append((f"fg:#ffffff bg:{ACCENT}", " " * max(52 - len(text), 1) + "\n"))
            else:
                colour = GOLD if kind in ("select", "type", "file") else ""
                out.append((f"fg:{colour}" if colour else "", f"     {text}\n"))
        if start + max_visible < total:
            out.append((f"fg:{DIM}", "     ↓ …\n"))
        return out

    def _metadata_lines(self):
        """Scrollable metadata block lines (header is drawn separately, fixed).

        Each line is either a single (style, text) fragment or a list of
        fragments (for parameter rows that mix colours).
        """
        lines = []
        for name, params in (self.metadata or []):
            lines.append((f"fg:{ACCENT} bold", f"   ▸ {name}"))
            for code, value in params:
                if not code:  # text block
                    for text_line in metadata_text_lines(value):
                        lines.append((f"fg:{DIM}", f"       {text_line}"))
                    continue
                desc = describe_param(code)
                row = [(f"fg:{GOLD}", f"       {code:<4}")]
                if desc:
                    row.append((f"fg:{DIM}", f"  {desc}"))
                row.append((f"fg:{DIM}", "  =  "))
                row.append(("", str(format_param_value(code, value))))
                lines.append(row)
        return lines

    def _body_metadata(self):
        if self.meta_error:
            return [
                ("", "\n"),
                (f"fg:{ACCENT} bold", f"   ✗ Could not read metadata\n\n"),
                (f"fg:{DIM}", f"   {self.meta_error}\n\n"),
                (f"fg:{DIM}", "   Press ← to pick another file.\n"),
            ]

        blocks = self.metadata or []
        nparams = sum(len(p) for _, p in blocks)

        # Fixed header (never scrolls): file name, counts, save confirmation.
        out = [("", "\n")]
        out.append((f"fg:{GOLD} bold", f"   {Path(self.selected_file).name}\n"))
        out.append((f"fg:{DIM}", f"   {len(blocks)} blocks · {nparams} parameters\n"))
        if self.meta_saved:
            out.append((f"fg:{GOLD} bold", f"   ✓ Saved to {self.meta_saved}\n"))
        else:
            out.append((f"fg:{DIM}", "   Press s to save all metadata to a .txt\n"))
        out.append(("", "\n"))

        # Scrollable block region.
        lines = self._metadata_lines()
        try:
            rows = self.app.output.get_size().rows
        except Exception:
            rows = 30
        max_visible = max(rows - 14, 6)  # room for banner, header and status
        total = len(lines)
        self.meta_scroll = max(0, min(self.meta_scroll, max(total - max_visible, 0)))
        start = self.meta_scroll

        out.append((f"fg:{DIM}", f"   ↑ more above ({start})\n") if start > 0 else ("", "\n"))
        for line in lines[start:start + max_visible]:
            if isinstance(line, list):
                out.extend(line)
                out.append(("", "\n"))
            else:
                style, text = line
                out.append((style, text + "\n"))
        if start + max_visible < total:
            out.append((f"fg:{DIM}", f"   ↓ more below ({total - start - max_visible})\n"))
        return out

    def _body_format(self):
        out = [("", "\n")]
        out.append((f"fg:{GOLD} bold", "   Which output formats would you like to generate?\n\n"))
        for i, (name, desc, _value) in enumerate(FORMAT_OPTIONS):
            selected = i == self.format_index
            if selected:
                out.append((f"fg:#ffffff bg:{ACCENT} bold", f"   › {name}"))
                out.append((f"fg:#ffffff bg:{ACCENT}", " " * max(28 - len(name), 1)))
                out.append((f"fg:#ffffff bg:{ACCENT}", f"{desc}"))
                out.append((f"fg:#ffffff bg:{ACCENT}", " " * 2 + "\n"))
            else:
                out.append(("", f"     {name}"))
                out.append(("", " " * max(28 - len(name), 1)))
                out.append((f"fg:{DIM}", f"{desc}\n"))
        return out

    def _body_confirm(self):
        out = [("", "\n")]
        formats = " and ".join(f".{f}" for f in self.selected_formats)
        out.append((f"fg:{DIM}", "   Folder    "))
        out.append(("", str(self.selected_folder) + "\n"))
        out.append((f"fg:{DIM}", "   Formats   "))
        out.append(("", formats + "\n\n"))

        if self.scanning:
            out.append((f"fg:{ACCENT} bold", f"   {self._spinner} "))
            out.append(("", "Scanning for OPUS files…\n"))
        elif not self.opus_files:
            out.append((f"fg:{ACCENT} bold", "   ✗  No OPUS files found in this folder.\n"))
            out.append((f"fg:{DIM}", "      Press ← to choose a different folder.\n"))
        else:
            out.append((f"fg:{GOLD} bold", f"   ✓  Found {len(self.opus_files):,} OPUS file(s).\n\n"))
            out.append((f"fg:{GOLD} bold", "   Press ⏎ to start converting.\n"))
        return out

    def _body_progress(self):
        out = [("", "\n")]
        total = len(self.opus_files)
        cur = self.done_count
        pct = (cur / total * 100) if total else 0
        width = 42
        filled = int(width * cur / total) if total else 0
        bar_done = "█" * filled
        bar_rest = "░" * (width - filled)

        eta = None
        if self.start_time and cur > 0 and cur < total:
            elapsed = time.time() - self.start_time
            eta = elapsed * total / cur - elapsed

        out.append((f"fg:{GOLD} bold", f"   {self._spinner} Converting…\n\n"))
        out.append(("   ", "   "))
        out.append((f"fg:{ACCENT}", bar_done))
        out.append((f"fg:{DIM}", bar_rest))
        out.append(("", f"  {pct:5.1f}%\n\n"))
        out.append((f"fg:{DIM}", "   files     "))
        out.append(("", f"{cur:,} / {total:,}\n"))
        out.append((f"fg:{DIM}", "   ETA       "))
        out.append(("", f"{_fmt_duration(eta)}\n"))
        if self.failed:
            out.append((f"fg:{ACCENT}", f"   failed    {self.failed:,}\n"))
        return out

    def _body_summary(self):
        out = [("", "\n")]
        total = self.success + self.failed
        formats = " and ".join(f".{f}" for f in self.selected_formats)
        out.append((f"fg:{GOLD} bold", "   Conversion summary\n\n"))
        out.append((f"fg:{DIM}", "   formats   "))
        out.append(("", formats + "\n"))
        out.append((f"fg:{DIM}", "   time      "))
        out.append(("", _fmt_duration(self.total_time) + "\n"))
        if self.success:
            avg = self.total_time / self.success if self.total_time else None
            out.append((f"fg:{DIM}", "   avg/file  "))
            out.append(("", _fmt_duration(avg) + "\n"))
        out.append(("", "\n"))
        out.append((f"fg:{GOLD} bold", f"   ✓  Converted   {self.success:,}\n"))
        if self.failed:
            rate = self.success / total * 100 if total else 0
            out.append((f"fg:{ACCENT} bold", f"   ✗  Failed      {self.failed:,}"))
            out.append((f"fg:{DIM}", f"   ({rate:.1f}% success)\n"))
        else:
            out.append((f"fg:{GOLD}", "   🎉  All files converted successfully!\n"))
        out.append(("", "\n"))
        out.append((f"fg:{DIM}", "   Press ⏎ or q to exit.\n"))
        return out

    def render_status(self):
        if self.step == "folder" and self._is_typing_path():
            hint = "type or paste a path   ⏎ go   esc cancel"
        elif self.step == "folder" and self._is_filtering():
            hint = "type to filter   ⏎ done   esc clear"
        elif self.step == "folder":
            hint = ("↑↓ move   ⏎ open   /  filter   ⇥ path   ← back   q quit"
                    if self.mode == "inspect"
                    else "↑↓ move   ⏎ open / select   ⇥ type path   ← back   q quit")
        elif self.step == "metadata":
            hint = ("← back   q quit" if self.meta_error
                    else "↑↓ scroll   s  save .txt   ← back   q quit")
        else:
            hint = {
                "welcome": "⏎ begin      q quit",
                "mode": "↑↓ move   ⏎ choose   ← back   q quit",
                "format": "↑↓ move   ⏎ choose   ← back   q quit",
                "confirm": "⏎ convert   ← back   q quit",
                "progress": "converting…  please wait",
                "summary": "⏎ / q  exit",
            }.get(self.step, "")
        return [("", "  OPUS dei  ·  "), ("bold", hint)]

    # -- Navigation --------------------------------------------------------

    def _move(self, delta):
        if self.step == "folder":
            n = len(self._visible_entries())
            if n:
                self.entry_index = (self.entry_index + delta) % n
        elif self.step == "format":
            self.format_index = (self.format_index + delta) % len(FORMAT_OPTIONS)
        elif self.step == "mode":
            self.mode_index = (self.mode_index + delta) % len(MODE_OPTIONS)
        elif self.step == "metadata":
            self.meta_scroll = max(0, self.meta_scroll + delta)

    def _enter_folder(self):
        """Reset the browser at the current directory for the active mode."""
        self.current_dir = Path.cwd()
        self.filter_field.text = ""
        self.path_error = None
        self._build_entries()
        self.step = "folder"

    def _activate(self):
        if self.step == "welcome":
            self.mode_index = 0
            self.step = "mode"

        elif self.step == "mode":
            self.mode = MODE_OPTIONS[self.mode_index][2]
            self._enter_folder()

        elif self.step == "folder":
            entries = self._visible_entries()
            if not entries:
                return
            self.entry_index = min(self.entry_index, len(entries) - 1)
            kind, _label, path = entries[self.entry_index]
            if kind == "select":
                self.selected_folder = str(self.current_dir)
                self.format_index = 0
                self.step = "format"
            elif kind == "type":
                self._focus_path_input()
            elif kind == "file":
                self._open_file(path)
            else:  # 'parent' or 'dir'
                self.current_dir = path
                self._build_entries()

        elif self.step == "format":
            self.selected_formats = FORMAT_OPTIONS[self.format_index][2]
            self.step = "confirm"
            self._start_scan()

        elif self.step == "confirm":
            if not self.scanning and self.opus_files:
                self.step = "progress"
                self._start_convert()

        elif self.step == "summary":
            self.app.exit()

    def _back(self):
        if self.step == "mode":
            self.step = "welcome"
        elif self.step == "folder":
            self.step = "mode"
        elif self.step == "format":
            self.step = "folder"
        elif self.step == "confirm":
            if not self.scanning:
                self.step = "format"
        elif self.step == "metadata":
            self.step = "folder"

    # -- Metadata inspector -------------------------------------------------

    def _open_file(self, path):
        """Parse an OPUS file's metadata and show the inspector."""
        self.selected_file = str(path)
        self.meta_scroll = 0
        self.meta_saved = None
        self.meta_error = None
        try:
            self.metadata = collect_metadata(str(path))
        except Exception as error:
            self.metadata = None
            self.meta_error = str(error) or error.__class__.__name__
        self.step = "metadata"

    def _save_metadata(self):
        """Write the current file's metadata to a .txt next to it."""
        if self.step != "metadata" or not self.metadata:
            return
        out_path = self.selected_file + ".metadata.txt"
        try:
            with open(out_path, "w", encoding="utf-8") as handle:
                handle.write(format_metadata_text(self.selected_file, self.metadata))
            self.meta_saved = out_path
        except Exception as error:
            self.meta_saved = f"(could not save: {error})"

    # -- Background workers -------------------------------------------------

    def _start_scan(self):
        self.scanning = True
        self.opus_files = []
        threading.Thread(target=self._scan_worker, daemon=True).start()

    def _scan_worker(self):
        try:
            files = find_opus_files(self.selected_folder)
        except Exception:
            files = []
        self.opus_files = files
        self.scanning = False

    def _start_convert(self):
        self.converting = True
        self.done_count = 0
        self.success = 0
        self.failed = 0
        self.start_time = time.time()
        threading.Thread(target=self._convert_worker, daemon=True).start()

    def _convert_worker(self):
        for i, opus_file in enumerate(self.opus_files):
            # Silence any parser warnings so they can't corrupt the TUI screen.
            with contextlib.redirect_stdout(io.StringIO()):
                try:
                    ok = convert_opus_file(opus_file, self.selected_formats)
                except Exception:
                    ok = False
            if ok:
                self.success += 1
            else:
                self.failed += 1
            self.done_count = i + 1
        self.total_time = time.time() - self.start_time
        self.converting = False
        self.step = "summary"

    # -- Application wiring -------------------------------------------------

    def _focus_filter_input(self):
        """Move focus to the live filter field (keeps existing filter text)."""
        self.app.layout.focus(self.filter_field)

    def _accept_filter(self, buffer):
        """Enter in the filter field just returns to the list (filter stays)."""
        self.app.layout.focus(self.body_window)
        return True  # keep the filter text applied

    def _build_application(self):
        # Direct path-entry field, shown during the folder step.
        self.path_field = TextArea(
            multiline=False,
            prompt=[("class:pathprompt", "   Path › ")],
            accept_handler=self._accept_path,
            style="class:pathinput",
            height=1,
        )
        # Live list-filter field, shown in the folder step of inspect mode.
        self.filter_field = TextArea(
            multiline=False,
            prompt=[("class:pathprompt", "   Filter › ")],
            accept_handler=self._accept_filter,
            style="class:pathinput",
            height=1,
        )
        # Reset the selection to the top whenever the filter text changes.
        self.filter_field.buffer.on_text_changed += lambda _: setattr(self, "entry_index", 0)

        typing = has_focus(self.path_field) | has_focus(self.filter_field)
        # While a text field is focused, letters (q, j, k, s, …) must type into
        # it, so navigation shortcuts are disabled. Arrow keys stay active so
        # the list can still be moved while filtering.
        nav = ~typing

        kb = KeyBindings()

        @kb.add("up")
        @kb.add("k", filter=nav)
        def _(event):
            self._move(-1)

        @kb.add("down")
        @kb.add("j", filter=nav)
        def _(event):
            self._move(1)

        @kb.add("enter", filter=nav)
        @kb.add("right", filter=nav)
        def _(event):
            self._activate()

        @kb.add("left", filter=nav)
        @kb.add("backspace", filter=nav)
        def _(event):
            self._back()

        @kb.add("q", filter=nav)
        def _(event):
            event.app.exit()

        @kb.add("c-c")
        def _(event):
            event.app.exit()

        # Save metadata to a .txt (metadata step only).
        @kb.add("s", filter=nav & Condition(lambda: self.step == "metadata"))
        def _(event):
            self._save_metadata()

        # Enter direct path-entry / filter modes (folder step, while browsing).
        folder_browsing = Condition(lambda: self.step == "folder") & nav

        @kb.add("tab", filter=folder_browsing)
        def _(event):
            self._focus_path_input()

        @kb.add("/", filter=folder_browsing & Condition(lambda: self.mode == "inspect"))
        def _(event):
            self._focus_filter_input()

        # Leave the path field without navigating.
        @kb.add("escape", filter=has_focus(self.path_field))
        def _(event):
            self.path_error = None
            event.app.layout.focus(self.body_window)

        # Leave the filter field and clear the filter.
        @kb.add("escape", filter=has_focus(self.filter_field))
        def _(event):
            self.filter_field.text = ""
            self.entry_index = 0
            event.app.layout.focus(self.body_window)

        self.body_window = Window(
            FormattedTextControl(self.render_body, focusable=True),
            dont_extend_height=True,
        )

        on_welcome = Condition(lambda: self.step == "welcome")
        on_folder = Condition(lambda: self.step == "folder")
        on_inspect_folder = Condition(lambda: self.step == "folder" and self.mode == "inspect")

        # Empty Window() fillers grow to absorb spare vertical space. The top
        # filler exists only on the welcome step (centring the splash); the
        # bottom filler is always present (pinning the status bar to the bottom
        # edge). Everything is a flat child of the root HSplit — a nested HSplit
        # would swallow the spare space instead of letting the fillers grow.
        root = HSplit([
            Window(height=self._welcome_top_height),  # exact top pad (welcome)
            Window(FormattedTextControl(self.render_title), height=len(TITLE_ART),
                   align=WindowAlign.CENTER),
            Window(FormattedTextControl(self.render_subtitle), height=1,
                   align=WindowAlign.CENTER),
            Window(height=1),
            ConditionalContainer(
                Window(FormattedTextControl(self.render_welcome_hint),
                       height=1, align=WindowAlign.CENTER),
                filter=on_welcome,
            ),
            self.body_window,
            ConditionalContainer(self.filter_field, filter=on_inspect_folder),
            ConditionalContainer(self.path_field, filter=on_folder),
            Window(),  # bottom filler
            Window(FormattedTextControl(self.render_status), height=1, style="class:status"),
        ])

        style = PTStyle.from_dict({
            "status": "reverse",
            "pathprompt": f"fg:{GOLD} bold",
            "pathinput": f"fg:#ffffff bg:#2a2a33",
        })

        return Application(
            layout=Layout(root),
            key_bindings=kb,
            style=style,
            full_screen=True,
            mouse_support=False,
            refresh_interval=0.1,  # keeps the spinner and progress bar animating
        )

    def run(self):
        self.app.run()


def main():
    """Run the OPUS dei converter as a full-screen TUI."""
    OpusDeiApp().run()


if __name__ == "__main__":
    main()