def process_j_new_blocks(filename):
    """
    Opens a file, reads it line by line, and processes blocks of lines
    starting with "j_new". For each block, it counts the occurrences
    of the number following "j_new" and prints the resulting list.
    The counts are reset after each block of "j_new" lines.

    Args:
        filename (str): The path to the input file.
    """
    try:
        with open(filename, 'r') as file:
            #j_new_counts = [0] * 8  # Initialize a list of 8 elements to zeros
            j_new_counts = [0] * 14  # Initialize a list of 14 elements to zeros
            in_j_new_block = False  # Flag to track if we are inside a j_new block

            for line_num, line in enumerate(file, 1):
                line_stripped = line.strip()

                if line_stripped.startswith("j_new"):
                    in_j_new_block = True
                    try:
                        # Attempt to extract the number after "j_new ="
                        # Assuming format like "j_new = X" or "j_new=X"
                        parts = line_stripped.split('=')
                        if len(parts) > 1:
                            number_str = parts[1].strip()
                            number = int(number_str)

                            # Increment the element corresponding to the number (1-based index)
                            #if 1 <= number <= 8:
                            if 0 <= number <= 13:  
                                j_new_counts[number] += 1
                            #else:
                                #print(f"Warning: Line {line_num}: Number '{number}' out of expected range (1-8). Skipping.")
                        else:
                            print(f"Warning: Line {line_num}: 'j_new' found but no number to extract. Skipping.")
                    except ValueError:
                        print(f"Warning: Line {line_num}: Could not parse number from '{line_stripped}'. Skipping.")
                else:
                    # If we were in a j_new block and now encounter a non-j_new line
                    if in_j_new_block:
                        print(f"\n--- End of j_new block (before line {line_num}) ---")
                        print("Counts for this block:")
                        for i, count in enumerate(j_new_counts):
                            print(f"  Index {i+1}: {count}")
                        
                        # Zero all values of the list for the next block
                        #j_new_counts = [0] * 8
                        j_new_counts = [0] * 14
                        in_j_new_block = False
                    # If not in a j_new block, just continue (ignore other lines)
        
            # After the loop, check if the file ended within a j_new block
            if in_j_new_block:
                print("\n--- End of j_new block (at EOF) ---")
                print("Counts for this block:")
                for i, count in enumerate(j_new_counts):
                    print(f"  Index {i+1}: {count}")

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def main():
    """
    Main function to get the filename from the user and start processing.
    """
    filename = input("Enter the name of the file to process: ")
    process_j_new_blocks(filename)

if __name__ == "__main__":
    main()