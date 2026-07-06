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
from pathlib import Path
from colorama import init, Fore, Back, Style
import time

from prompt_toolkit.application import Application
from prompt_toolkit.key_binding import KeyBindings
from prompt_toolkit.layout import Layout
from prompt_toolkit.layout.containers import HSplit, Window, ConditionalContainer
from prompt_toolkit.layout.controls import FormattedTextControl
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

# Little mascot kept from previous versions — the OPUS dei signature.
MASCOT_ART = r"""
                     ___
               __   /   \  _____
              /(o)\/     \/    /\
              \___/\     /\___/  \
             /     /\___/     \ / \
            /     / / /\ \     \ / \
           /     / / /  \ \     \  /
          /     ^ ^ ^    ^ ^     \/
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
        wavenumber_step = -(first_wavenumber - last_wavenumber) / len(absorption_spectrum)

        # Create the wavenumber axis
        wavenumbers = np.arange(first_wavenumber, last_wavenumber, wavenumber_step)

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


# Output format options presented in the format-selection step.
FORMAT_OPTIONS = [
    ("Both .dpt and .mzz", "full resolution + compressed archive", ["dpt", "mzz"]),
    ("Only .dpt", "full resolution tab-delimited text", ["dpt"]),
    ("Only .mzz", "compressed 1 cm⁻¹ ZIP archive", ["mzz"]),
]

SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"


class OpusDeiApp:
    """A full-screen, navigable TUI for the OPUS dei converter.

    The interface is a fixed screen (banner on top, a body that changes with
    the current step, and a status bar with shortcuts at the bottom). Nothing
    scrolls — the whole screen repaints in place, giving a smooth, modern feel.
    The spectral parsing/conversion engine is untouched and runs in a
    background thread so the UI stays responsive.
    """

    def __init__(self):
        self.step = "welcome"          # welcome | folder | format | confirm | progress | summary

        # Folder browser state.
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
        """List selectable items for the folder browser at the current path."""
        entries = [("select", f"✓  Use this folder", self.current_dir)]
        parent = self.current_dir.parent
        if parent != self.current_dir:
            entries.append(("parent", "..  (parent folder)", parent))
        try:
            subdirs = sorted(
                (p for p in self.current_dir.iterdir() if p.is_dir()),
                key=lambda p: p.name.lower(),
            )
        except (PermissionError, OSError):
            subdirs = []
        for p in subdirs:
            entries.append(("dir", p.name + "/", p))
        self.entries = entries
        self.entry_index = 0

    # -- Rendering ---------------------------------------------------------

    def render_title(self):
        frags = title_fragments()
        # Drop the final newline so the art fits an exact-height window.
        style, text = frags[-1]
        frags[-1] = (style, text.rstrip("\n"))
        return frags

    def render_subtitle(self):
        return [
            (f"fg:{GOLD}", "   Bruker OPUS spectral converter"),
            (f"fg:{DIM}", "   ·   "),
            (f"fg:{GOLD}", "v3.0"),
            (f"fg:{DIM}", "   ·   Mario González-Jiménez · University of Glasgow"),
        ]

    def render_body(self):
        return getattr(self, f"_body_{self.step}")()

    def _body_welcome(self):
        out = [("", "\n")]
        rows = [
            ("What it does", ".dpt / .mzz from Bruker OPUS files"),
            ("Recursive", "finds files in every subfolder"),
            ("Batch", "handles dozens to tens of thousands"),
        ]
        inner = 52
        out.append((f"fg:{ACCENT}", "   ╭" + "─" * inner + "╮\n"))
        for label, value in rows:
            line = f" {label}   {value}"
            pad = " " * max(inner - len(line), 0)
            out.append((f"fg:{ACCENT}", "   │"))
            out.append((f"fg:{DIM}", f" {label}   "))
            out.append(("", value))
            out.append(("", pad))
            out.append((f"fg:{ACCENT}", "│\n"))
        out.append((f"fg:{ACCENT}", "   ╰" + "─" * inner + "╯\n"))
        out.append(("", "\n"))
        for line in MASCOT_ART:
            out.append((f"fg:{DIM}", "     " + line + "\n"))
        out.append(("", "\n"))
        out.append((f"fg:{GOLD} bold", "   Press ⏎ to begin.\n"))
        return out

    def _is_typing_path(self):
        """True while the direct path-entry field holds keyboard focus."""
        try:
            return self.app.layout.has_focus(self.path_field)
        except Exception:
            return False

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
        if resolved.is_dir():
            target = resolved
        elif resolved.is_file():
            target = resolved.parent  # accept a file path, use its folder
        else:
            self.path_error = f"Not a folder: {raw}"
            return True  # keep the text so the user can fix it
        self.current_dir = target
        self._build_entries()
        self.path_error = None
        self.app.layout.focus(self.body_window)  # back to the browser
        return False

    def _body_folder(self):
        out = [("", "\n")]
        out.append((f"fg:{GOLD} bold", "   Select the folder containing your OPUS files\n\n"))
        out.append((f"fg:{DIM}", "   " + str(self.current_dir) + "\n"))
        if self.path_error:
            out.append((f"fg:{ACCENT}", "   ✗ " + self.path_error + "\n"))
        out.append(("", "\n"))

        # Simple viewport so long directory listings never overflow.
        max_visible = 12
        total = len(self.entries)
        start = max(0, min(self.entry_index - max_visible // 2, total - max_visible))
        start = max(start, 0)
        visible = self.entries[start:start + max_visible]

        if start > 0:
            out.append((f"fg:{DIM}", "     ↑ …\n"))
        for i, (kind, label, _path) in enumerate(visible, start=start):
            selected = i == self.entry_index
            icon = {"select": "✓ ", "parent": "↩ ", "dir": "📁 "}[kind]
            if kind == "select":
                icon = "✓ "
            text = f"{icon}{label}"
            if selected:
                out.append((f"fg:#ffffff bg:{ACCENT} bold", f"   › {text}"))
                out.append((f"fg:#ffffff bg:{ACCENT}", " " * max(52 - len(text), 1) + "\n"))
            else:
                colour = GOLD if kind == "select" else ""
                out.append((f"fg:{colour}" if colour else "", f"     {text}\n"))
        if start + max_visible < total:
            out.append((f"fg:{DIM}", "     ↓ …\n"))
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
        else:
            hints = {
                "welcome": "⏎ begin      q quit",
                "folder": "↑↓ move   ⏎ open / select   ⇥ type path   ← back   q quit",
                "format": "↑↓ move   ⏎ choose   ← back   q quit",
                "confirm": "⏎ convert   ← back   q quit",
                "progress": "converting…  please wait",
                "summary": "⏎ / q  exit",
            }
            hint = hints.get(self.step, "")
        return [("", "  OPUS dei  ·  "), ("bold", hint)]

    # -- Navigation --------------------------------------------------------

    def _move(self, delta):
        if self.step == "folder" and self.entries:
            self.entry_index = (self.entry_index + delta) % len(self.entries)
        elif self.step == "format":
            self.format_index = (self.format_index + delta) % len(FORMAT_OPTIONS)

    def _activate(self):
        if self.step == "welcome":
            self.current_dir = Path.cwd()
            self._build_entries()
            self.step = "folder"

        elif self.step == "folder":
            if not self.entries:
                return
            kind, _label, path = self.entries[self.entry_index]
            if kind == "select":
                self.selected_folder = str(self.current_dir)
                self.format_index = 0
                self.step = "format"
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
        if self.step == "folder":
            self.step = "welcome"
        elif self.step == "format":
            self.step = "folder"
        elif self.step == "confirm":
            if not self.scanning:
                self.step = "format"

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

    def _build_application(self):
        # Direct path-entry field, shown only during the folder step.
        self.path_field = TextArea(
            multiline=False,
            prompt=[("class:pathprompt", "   path › ")],
            accept_handler=self._accept_path,
            style="class:pathinput",
            height=1,
        )
        typing = has_focus(self.path_field)
        # Navigation keys are inactive while the path field is focused, so
        # letters (q, j, k, /, …) type into the path instead of triggering nav.
        nav = ~typing

        kb = KeyBindings()

        @kb.add("up", filter=nav)
        @kb.add("k", filter=nav)
        def _(event):
            self._move(-1)

        @kb.add("down", filter=nav)
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

        # Enter direct path-entry mode (folder step only, when browsing).
        folder_browsing = Condition(lambda: self.step == "folder") & nav

        @kb.add("tab", filter=folder_browsing)
        @kb.add("/", filter=folder_browsing)
        def _(event):
            self.path_error = None
            self.path_field.text = ""  # start empty so a pasted path just works
            event.app.layout.focus(self.path_field)

        # Leave path-entry mode without navigating.
        @kb.add("escape", filter=typing)
        def _(event):
            self.path_error = None
            event.app.layout.focus(self.body_window)

        self.body_window = Window(FormattedTextControl(self.render_body, focusable=True))

        root = HSplit([
            Window(height=1),
            Window(FormattedTextControl(self.render_title), height=len(TITLE_ART)),
            Window(FormattedTextControl(self.render_subtitle), height=1),
            Window(height=1),
            self.body_window,
            ConditionalContainer(
                self.path_field,
                filter=Condition(lambda: self.step == "folder"),
            ),
            Window(FormattedTextControl(self.render_status), height=1, style="class:status"),
        ])

        style = PTStyle.from_dict({
            "status": "reverse",
            "pathprompt": f"fg:{GOLD} bold",
            "pathinput": f"fg:#ffffff bg:#2a2a33",
        })

        return Application(
            layout=Layout(root, focused_element=self.body_window),
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