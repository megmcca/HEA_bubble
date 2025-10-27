import io

def read_file_content(filename):
    """
    Reads the cluster file and returns the number of site of spins from 1 to 8.

    Args:
        filename (str): The name of the file to read.

    Returns:
        for each time, the number of sites with that spin.
        
    """
    all_processed_blocks = []
    try:
        with open(filename, "r") as file:
            while True:
                current_time_value = None
                header_found_for_block = False
                spin_number_list = [0.0] * 8
                data_lines_in_block = False # Flag to check if any data was processed in this block
                size_3_to_8_total = 0 # Variable to store the sum of sizes for spins 3 to 8

                # --- Phase 1: Find "Time =" line and header for the current block ---
                # Read lines until "Time =" or "id ..." or EOF
                line = file.readline()
                while line:
                    line_stripped = line.strip()
                    if "Time = " in line_stripped:
                        try:
                            current_time_value = float(line_stripped.split("=")[1].strip())
                        except ValueError:
                            print(f"Warning: Could not extract time value from line: {line_stripped}")
                            current_time_value = None
                    elif "id ivalue dvalue size cx cy cz xlo xhi ylo yhi zlo zhi" in line_stripped:
                        header_found_for_block = True
                        break # Found header, proceed to data processing
                    
                    line = file.readline() # Read next line

                if not line: # EOF reached while looking for Time or Header
                    break # Exit the main processing loop

                if not header_found_for_block:
                    # This case means "Time =" was found, but no header followed before EOF or next block
                    # Or, if Time was not found, but we hit EOF
                    break # No valid block start, exit

                # --- Phase 2: Process data lines for the current block ---
                line = file.readline() # Read the first data line after the header
                while line:
                    line_stripped = line.strip()
                    if not line_stripped: # Blank line indicates end of data block
                        break
                    
                    try:
                        values = line_stripped.split()
                        if len(values) >= 4:
                            data_lines_in_block = True
                            ivalue = int(values[1])
                            size = int(values[3])

                            # Validate ivalue for list indexing (must be 1 to 8)
                            if 1 <= ivalue <= 8:
                                spin_number_list[ivalue - 1] += size
                                if 3 <= ivalue <= 8:
                                    size_3_to_8_total += size   
                            else:
                                print(f"Warning: ivalue {ivalue} out of expected range (1-8) for line: {line_stripped}")
                        else:
                            print(f"Warning: Line has fewer than 4 values, skipping: {line_stripped}")
                    except ValueError:
                        print(f"Warning: Could not parse numbers from line, skipping: {line_stripped}")
                    except IndexError:
                        print(f"Error: Indexing issue with ivalue for line: {line_stripped}")
                    
                    line = file.readline() # Read next line for this block

                # --- Phase 3: Store processed block data ---
                # Only store if a time was found and data was actually processed in this block
                if current_time_value is not None and data_lines_in_block:
                    all_processed_blocks.append((current_time_value, spin_number_list,  size_3_to_8_total))
                
                # If the inner loop broke because of EOF, then the outer loop should also break
                if not line and line_stripped: # EOF reached, and it was not a blank line
                    break
                elif not line and not line_stripped: # EOF reached, and it was a blank line
                    break


        return all_processed_blocks

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return []
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        return []

def main():
    """
    Main function to get the filename from the user, read and process all data blocks,
    and print the extracted time and calculated lists for each block.
    """
    filename = input("Enter the name of the file to read: ")
    output_filename = f"{filename}.size" # Create the output filename
    all_data_blocks = read_file_content(filename)
    
    if all_data_blocks:
        print("\n--- Processed Data Blocks ---")
        try:
            with open(output_filename, "w") as outfile:
                outfile.write("--- Processed Data Blocks ---\n")
                for i, (time_value, Fspin_number_list, size_3_to_8_total) in enumerate(all_data_blocks):
                    print(f"Time = {time_value}")
                    print(f"Spin    Size")
                    for j, total_size in enumerate(Fspin_number_list):
                        print(f"{j+1}, {total_size}")
                    print(f"Size of spins 3 to 8: {size_3_to_8_total}") # Print the
                    print("\n")

                    # Write to file
                    outfile.write(f"Time = {time_value}\n")
                    outfile.write(f"Spin     Size\n")
                    for j, total_size in enumerate(Fspin_number_list):
                        outfile.write(f"{j+1}, {total_size}\n")
                    outfile.write(f"Size of spins 3 to 8: {size_3_to_8_total}\n")
                    outfile.write("\n")
                print(f"\nOutput also saved to '{output_filename}'")
        except Exception as e:
            print(f"Error writing to output file '{output_filename}': {e}")    
    else:
        print("No data blocks were processed or an error occurred.")

if __name__ == "__main__":
    main()

