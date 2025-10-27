import os

def read_and_process_file(input_filename, output_filename):
    try:
        with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
            reading_atoms = False
            atom_id_counter = 0  # Initialize the counter
            for line in input_file:
                line = line.strip()
                if "ITEM: TIMESTEP" in line:
                    output_file.write("ITEM: TIMESTEP\n")
                    next_line = input_file.readline().strip()
                    output_file.write(f"{next_line}\n")
                    reading_atoms = False
                elif "ITEM: NUMBER OF ATOMS" in line:
                    output_file.write("ITEM: NUMBER OF ATOMS\n")
                    number_of_atoms_line = input_file.readline().strip()
                    output_file.write(f"{number_of_atoms_line}\n")
                    reading_atoms = False
                elif "ITEM: BOX BOUNDS" in line:
                    output_file.write("ITEM: BOX BOUNDS\n")
                    for _ in range(3):
                        box_bounds_line = input_file.readline().strip()
                        output_file.write(f"{box_bounds_line}\n")
                    reading_atoms = False
                    output_file.write("ITEM: ATOMS id i1 x y z\n")
                elif "ITEM: ATOMS" in line:
                    reading_atoms = True
                    atom_id_counter = 1 # Reset counter at the start of new ITEM: ATOMS section
                elif reading_atoms:
                    numbers_str = line.split()
                    if numbers_str:
                        try:
                            numbers = [float(num) for num in numbers_str]
                            if len(numbers) >= 2:
                                if numbers[1] == 3:
                                    if len(numbers) >= 4:
                                        modified_numbers_add = list(numbers)
                                        modified_numbers_subtract = list(numbers)

                                        modified_numbers_add[2] += 0.2
                                        modified_numbers_add[3] += 0.2
                                        modified_numbers_subtract[2] -= 0.2
                                        modified_numbers_subtract[3] -= 0.2
                                        
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_add[2]:.6f} {modified_numbers_add[3]:.6f} {' '.join(map(str, modified_numbers_add[4:]))}\n")
                                        atom_id_counter += 1
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_subtract[2]:.6f} {modified_numbers_subtract[3]:.6f} {' '.join(map(str, modified_numbers_subtract[4:]))}\n")
                                        atom_id_counter += 1
                                elif numbers[1] == 4:
                                    if len(numbers) >= 5:
                                        modified_numbers_add = list(numbers)
                                        modified_numbers_subtract = list(numbers)
                                        modified_numbers_add[3] += 0.2
                                        modified_numbers_add[4] += 0.2
                                        modified_numbers_subtract[3] -= 0.2
                                        modified_numbers_subtract[4] -= 0.2
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_add[2]:.6f} {modified_numbers_add[3]:.6f} {modified_numbers_add[4]:.6f} {' '.join(map(str, modified_numbers_add[5:]))}\n")
                                        atom_id_counter += 1
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_subtract[2]:.6f} {modified_numbers_subtract[3]:.6f} {modified_numbers_subtract[4]:.6f} {' '.join(map(str, modified_numbers_subtract[5:]))}\n")
                                        atom_id_counter += 1
                                elif numbers[1] == 5:
                                    if len(numbers) >= 5:
                                        modified_numbers_add = list(numbers)
                                        modified_numbers_subtract = list(numbers)
                                        modified_numbers_add[2] += 0.2
                                        modified_numbers_add[4] += 0.2
                                        modified_numbers_subtract[2] -= 0.2
                                        modified_numbers_subtract[4] -= 0.2
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_add[2]:.6f} {modified_numbers_add[3]:.6f} {modified_numbers_add[4]:.6f} {' '.join(map(str, modified_numbers_add[5:]))}\n")
                                        atom_id_counter += 1
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_subtract[2]:.6f} {modified_numbers_subtract[3]:.6f} {modified_numbers_subtract[4]:.6f} {' '.join(map(str, modified_numbers_subtract[5:]))}\n")
                                        atom_id_counter += 1
                                elif numbers[1] == 6:
                                    if len(numbers) >= 5:
                                        modified_numbers_add = list(numbers)
                                        modified_numbers_subtract = list(numbers)
                                        modified_numbers_add[2] += 0.2
                                        modified_numbers_subtract[3] -= 0.2
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_add[2]:.6f} {modified_numbers_add[3]:.6f} {modified_numbers_add[4]:.6f} {' '.join(map(str, modified_numbers_add[5:]))}\n")
                                        atom_id_counter += 1
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_subtract[2]:.6f} {modified_numbers_subtract[3]:.6f} {modified_numbers_subtract[4]:.6f} {' '.join(map(str, modified_numbers_subtract[5:]))}\n")
                                        atom_id_counter += 1
                                elif numbers[1] == 7:
                                    if len(numbers) >= 5:
                                        modified_numbers_add = list(numbers)
                                        modified_numbers_subtract = list(numbers)
                                        modified_numbers_add[3] += 0.2
                                        modified_numbers_subtract[4] -= 0.2
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_add[2]:.6f} {modified_numbers_add[3]:.6f} {modified_numbers_add[4]:.6f} {' '.join(map(str, modified_numbers_add[5:]))}\n")
                                        atom_id_counter += 1
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_subtract[2]:.6f} {modified_numbers_subtract[3]:.6f} {modified_numbers_subtract[4]:.6f} {' '.join(map(str, modified_numbers_subtract[5:]))}\n")
                                        atom_id_counter += 1
                                elif numbers[1] == 8:
                                    if len(numbers) >= 5:
                                        modified_numbers_add = list(numbers)
                                        modified_numbers_subtract = list(numbers)
                                        modified_numbers_add[2] += 0.2
                                        modified_numbers_subtract[4] -= 0.2
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_add[2]:.6f} {modified_numbers_add[3]:.6f} {modified_numbers_add[4]:.6f} {' '.join(map(str, modified_numbers_add[5:]))}\n")
                                        atom_id_counter += 1
                                        output_file.write(f"{atom_id_counter} {int(numbers[1])} {modified_numbers_subtract[2]:.6f} {modified_numbers_subtract[3]:.6f} {modified_numbers_subtract[4]:.6f} {' '.join(map(str, modified_numbers_subtract[5:]))}\n")
                                        atom_id_counter += 1
                                else:
                                    output_file.write(f"{atom_id_counter} {int(numbers[1])} {' '.join(map(str, numbers[2:]))}\n")
                                    atom_id_counter += 1
                        except ValueError:
                            pass

    except FileNotFoundError:
        print(f"The file '{input_filename}' does not exist.")
    except PermissionError:
        print(f"Permission denied when trying to access '{input_filename}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    input_filename = os.path.expanduser('~/Documents/myLocalRepo/spparks_bcc_dumbbell-main/dump.diff')
    output_filename = 'dump.diff.exp'
    read_and_process_file(input_filename, output_filename)

if __name__ == "__main__":
    main()
