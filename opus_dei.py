#!/usr/bin/env python3
"""
OPUS File Converter

This script converts Bruker OPUS spectral files to two formats:
1. .dpt files: Tab-delimited text files with wavenumber and absorption data
2. .mzz files: Compressed format with rounded wavenumbers for space efficiency

Author: Converted from Jupyter notebook
Usage: python opus_converter.py
"""

import numpy as np
import struct
import os
import zipfile
import questionary
from pathlib import Path
from colorama import init, Fore, Back, Style
import time

# Initialise colorama for cross-platform colour support
init(autoreset=True)


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


def welcome_message():
    """
    Prints the welcome message for OPUS File Converter.
    """
    print(f"\n\n{Fore.YELLOW}**********************************************************************{Style.RESET_ALL}")
    print(f"{Fore.YELLOW}* Welcome to OPUS Converter — a Bruker OPUS file processing tool.   *{Style.RESET_ALL}")
    print(f"{Fore.YELLOW}**********************************************************************{Style.RESET_ALL}")
    
    banner = rf"""

                     ___
                    /   \         ____
              -----/ (o) \-------/    \
                   \___  /       \__  /              {Fore.LIGHTRED_EX}  OPUS dei v3.0{Style.RESET_ALL}
                       ||           ||               {Fore.LIGHTRED_EX}  Mario González-Jiménez{Style.RESET_ALL}
                      //\\          //\\             {Fore.LIGHTRED_EX}  22 Sept 2025 - University of Glasgow{Style.RESET_ALL}
                     //  \\        //  \\
                    ^^    ^^      ^^    ^^
            """
    print(banner)
    
    # Wait 3 seconds before showing the rest
    time.sleep(3)
    
    print(f"\n\n{Fore.YELLOW}The conversion process follows these steps:{Style.RESET_ALL}\n")
    
    print(f"  {Fore.GREEN}1. Directory Selection{Style.RESET_ALL}")
    print(f"     Browse and select the folder containing your OPUS files.")
    print(f"     The program will recursively search through all subdirectories")
    print(f"     to find files with numerical endings (e.g., spectrum.0, data.1).")
    print(f"     {Fore.CYAN}Interactive folder browser with full path support.{Style.RESET_ALL}\n")
    
    print(f"  {Fore.GREEN}2. Output Format Selection{Style.RESET_ALL}")
    print(f"     Choose your preferred output formats:")
    print(f"     • {Fore.BLUE}.dpt files{Style.RESET_ALL} - Full resolution tab-delimited text files")
    print(f"       Perfect for direct import into analysis software")
    print(f"     • {Fore.YELLOW}.mzz files{Style.RESET_ALL} - Compressed format with 1 cm⁻¹ resolution")
    print(f"       Space-efficient ZIP archives for long-term storage")
    print(f"     {Fore.CYAN}Mix and match formats based on your workflow needs.{Style.RESET_ALL}\n")
    
    print(f"  {Fore.GREEN}3. Batch Processing{Style.RESET_ALL}")
    print(f"     Intelligent processing with adaptive progress display:")
    print(f"     • Small batches (≤100 files): Individual file tracking")
    print(f"     • Large batches (>100 files): Progress bar with ETA")
    print(f"     • Error handling and detailed conversion statistics")
    print(f"     {Fore.CYAN}Optimised for datasets from dozens to tens of thousands of files.{Style.RESET_ALL}\n")
    
    print(f"\n{Fore.YELLOW}Ready to convert your OPUS files? Let's get started!{Style.RESET_ALL}")
    
    # Small pause before continuing
    time.sleep(2)


def print_header():
    """
    Print a colourful header for the application.
    """
    print(f"\n{Fore.MAGENTA}{'='*60}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{Style.BRIGHT}🔄 FILE CONVERSION INTERFACE{Style.RESET_ALL}")
    print(f"{Fore.MAGENTA}{'='*60}{Style.RESET_ALL}\n")


def print_progress_bar(current, total, width=50, start_time=None):
    """
    Print a progress bar with percentage, visual bar, and estimated time remaining.
    
    Args:
        current (int): Current progress count
        total (int): Total items to process
        width (int): Width of the progress bar in characters
        start_time (float): Time when processing started (for ETA calculation)
    """
    if total == 0:
        return
        
    percentage = (current / total) * 100
    filled_width = int(width * current // total)
    
    # Create the visual bar
    bar = f"{'█' * filled_width}{'░' * (width - filled_width)}"
    
    # Calculate ETA if start_time is provided
    eta_text = ""
    if start_time and current > 0:
        elapsed = time.time() - start_time
        if current < total:
            estimated_total = elapsed * total / current
            eta_seconds = estimated_total - elapsed
            if eta_seconds > 60:
                eta_text = f" | ETA: {int(eta_seconds // 60)}m {int(eta_seconds % 60)}s"
            else:
                eta_text = f" | ETA: {int(eta_seconds)}s"
        else:
            eta_text = f" | Completed in: {int(elapsed)}s"
    
    # Print the progress bar (using \r to overwrite the same line)
    print(f"\r{Fore.CYAN}Progress:{Style.RESET_ALL} [{bar}] {percentage:6.2f}% ({current:,}/{total:,}){eta_text}", end='', flush=True)


def print_summary(successful, failed, formats, total_time=None):
    """
    Print a colourful summary of the conversion results.
    
    Args:
        successful (int): Number of successful conversions
        failed (int): Number of failed conversions
        formats (list): List of output formats used
        total_time (float): Total time taken for processing
    """
    print(f"\n\n{Fore.MAGENTA}{'='*60}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{Style.BRIGHT}📊 CONVERSION SUMMARY{Style.RESET_ALL}")
    print(f"{Fore.MAGENTA}{'='*60}{Style.RESET_ALL}")
    
    format_names = " and ".join([f".{fmt}" for fmt in formats])
    print(f"{Fore.WHITE}Output format(s):{Style.RESET_ALL} {format_names}")
    
    if total_time:
        if total_time > 60:
            time_text = f"{int(total_time // 60)}m {int(total_time % 60)}s"
        else:
            time_text = f"{total_time:.1f}s"
        print(f"{Fore.WHITE}Total processing time:{Style.RESET_ALL} {time_text}")
        
        if successful > 0:
            avg_time = total_time / successful
            if avg_time < 1:
                avg_text = f"{avg_time*1000:.0f}ms"
            else:
                avg_text = f"{avg_time:.2f}s"
            print(f"{Fore.WHITE}Average time per file:{Style.RESET_ALL} {avg_text}")
    
    if successful > 0:
        print(f"{Fore.GREEN}✅ Successfully converted:{Style.RESET_ALL} {Style.BRIGHT}{successful:,}{Style.RESET_ALL} files")
    
    if failed > 0:
        print(f"{Fore.RED}❌ Failed conversions:{Style.RESET_ALL} {Style.BRIGHT}{failed:,}{Style.RESET_ALL} files")
        success_rate = (successful / (successful + failed)) * 100
        print(f"{Fore.YELLOW}📈 Success rate:{Style.RESET_ALL} {success_rate:.1f}%")
    else:
        print(f"{Fore.GREEN}🎉 All files converted successfully!{Style.RESET_ALL}")
    
    print(f"{Fore.MAGENTA}{'='*60}{Style.RESET_ALL}\n")


def main():
    """
    Main function to run the OPUS file conversion process.
    """
    # Show welcome message first
    welcome_message()
    
    print_header()
    
    # Get directory path from user
    folder_path = questionary.path(
        f"Please select the folder containing OPUS files:",
        only_directories=True,
        qmark="📂",
        style=questionary.Style([
            ('question', 'fg:#ff0066 bold'),
            ('answer', 'fg:#44ff00 bold'),
        ])
    ).ask()
    
    if not folder_path:
        print(f"{Fore.RED}❌ No folder selected. Exiting.{Style.RESET_ALL}")
        return
    
    # Ask user which output formats they want
    output_format_choice = questionary.select(
        f"🎯 Which output formats would you like to generate?",
        choices=[
            {"name": f"Both .dpt and .mzz files", "value": ["dpt", "mzz"]},
            {"name": f"Only .dpt files (full resolution)", "value": ["dpt"]},
            {"name": f"Only .mzz files (compressed)", "value": ["mzz"]}
        ],
        style=questionary.Style([
            ('question', 'fg:#ff0066 bold'),
            ('answer', 'fg:#44ff00 bold'),
            ('pointer', 'fg:#ff0066 bold'),
            ('highlighted', 'fg:#ff0066 bold'),
        ])
    ).ask()
    
    if not output_format_choice:
        print(f"{Fore.RED}❌ No format selected. Exiting.{Style.RESET_ALL}")
        return
    
    # Find all OPUS files
    print(f"{Fore.CYAN}🔍 Searching for OPUS files...{Style.RESET_ALL}")
    opus_files = find_opus_files(folder_path)
    
    if not opus_files:
        print(f"{Fore.RED}❌ No OPUS files found in {folder_path}{Style.RESET_ALL}")
        return
    
    format_names = " and ".join([f".{fmt}" for fmt in output_format_choice])
    print(f"{Fore.GREEN}✅ Found {Style.BRIGHT}{len(opus_files):,}{Style.RESET_ALL}{Fore.GREEN} OPUS file(s) to convert to {format_names} format(s).{Style.RESET_ALL}")
    
    # Determine if we should show individual file processing
    show_individual = len(opus_files) <= 100
    
    if not show_individual:
        print(f"{Fore.YELLOW}📋 Large batch detected - using progress bar mode{Style.RESET_ALL}")
    
    print(f"\n{Fore.MAGENTA}{'─'*60}{Style.RESET_ALL}")
    print(f"{Fore.CYAN}{Style.BRIGHT}🚀 STARTING CONVERSION{Style.RESET_ALL}")
    print(f"{Fore.MAGENTA}{'─'*60}{Style.RESET_ALL}")
    
    # Convert each file
    successful_conversions = 0
    failed_conversions = 0
    start_time = time.time()
    
    for i, opus_file in enumerate(opus_files):
        if not show_individual:
            # Show progress bar for large batches
            print_progress_bar(i, len(opus_files), start_time=start_time)
        else:
            # Show individual file processing for small batches
            print(f"{Fore.WHITE}[{i+1}/{len(opus_files)}]{Style.RESET_ALL}", end=" ")
            
        if convert_opus_file(opus_file, output_format_choice, show_individual_files=show_individual):
            successful_conversions += 1
        else:
            failed_conversions += 1
    
    # Final progress bar update
    if not show_individual:
        print_progress_bar(len(opus_files), len(opus_files), start_time=start_time)
    
    total_time = time.time() - start_time
    print_summary(successful_conversions, failed_conversions, output_format_choice, total_time)


if __name__ == "__main__":
    main()