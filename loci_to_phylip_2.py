#!/usr/local/bin/python3
import itertools
import sys
import numpy as np

def initialize_variables(args):
    """
    Orchestrates initialization of the global variables that will be used in this script.

    Args:
        args: The command line arguments passed to the program.
    """
    set_file_names(args)
    set_loci_groups()
    initialize_sequences()
    initialize_sequence_matrix()

def set_file_names(command_line_args):
    """
    Sets the names of the input and output files that will be used.

    Args:
        command_line_args: The command line arguments passed to the program.
    """
    global input_file_name
    global output_file_name
    global arguments_file_name
    global matrix_file_name
    if len(command_line_args) != 5:
        print("Missing one or more required arguments. Required arguments are:\n1. Input File name,\n" +
              "2. Arguments File name,\n3. Phylip Output File name,\n4. Matrix Output File name.")
        sys.exit(1)
    input_file_name = command_line_args[1]
    arguments_file_name = command_line_args[2]
    output_file_name = command_line_args[3]
    matrix_file_name = command_line_args[4]

def set_loci_groups():
    """
    Reads in the loci groups that should be extracted from the .loci file.
    """
    global loci_groups
    loci_groups = []
    reading_groups = False # we haven't reached the loci groups until we reach the word 'extract'
    for line in open(arguments_file_name, 'r'):
        if 'extract' in line.lower():
            reading_groups = True # now we start reading the groups
            continue
        elif not reading_groups:
            continue # skip lines until we reach the loci groups
        elif 'minimum of ' in line: # this is a line corresponding to a group
            groups_in_line = process_argument_line(line)
            loci_groups.extend(groups_in_line)
        else: # we've finished reading the groups
            break

def process_argument_line(line):
    """
    Processes a line of the "Extract" section of the argument file.

    Args:
        line: A line containing information about a loci group to be extracted.

    Returns:
        An array of loci groups of the format
        {'minimum_matches': 2, 'species': ('C10', 'F10', 'G10', 'H10', 'A11')}
    """
    groups_in_line = []
    min_and_species = line[11:].split('#')[0].rstrip() # remove 'minimum of ', comments, and trailing whitespace
    min = min_and_species[0]
    if not min.isdigit():
        raise ValueError('There was a problem with the following line of the arguments file:\n' + line +
                         '\nMake sure that loci extraction parameters follow the format "minimum of x (...)"' +
                         ' and that x is an integer.')
    species_string = min_and_species[3:-1] # leave out min and surrounding parentheses
    species = parse_species_object(species_string)
    for spec in species:
        # create a separate group for each permutation of species
        groups_in_line.append({'minimum_matches': int(min), 'species': spec})
    return groups_in_line

def parse_species_object(species_string):
    """
    Parses groups of species from a given collection string (i.e. "(A11 and/or B11, C2)")

    Args:
        species_string: A string representing a group of species.
    Returns:
        An array of tuples corresponding to combinations of given species.
    """
    requirements = species_string.split(', ')
    nested_requirements = [x.split(' and/or ') for x in requirements]
    if len(nested_requirements) == 1: # species list has no commas, just and/or's
        return [tuple(nested_requirements[0])] # so we don't need permutations
    else:
        return list(itertools.product(*nested_requirements))

def initialize_sequences():
    """
    Initializes the hash of id-sequence pairs that will keep track of all of the sequences
    as they are being built.
    Keys are initialized here, values are empty strings.
    """
    global sequences
    sequences = {}
    reading_ids = False # we haven't reached the ids until we reach a line containing the word 'id'
    for line in open(arguments_file_name, 'r'):
        if 'id' in line.lower():
            reading_ids = True
            continue
        elif not reading_ids:
            continue
        elif len(line) >= 2:
            sequences[line.rstrip()] = ''
        else:
            break

def initialize_sequence_matrix():
    """
    Initializes the matrix of the number of times species occur in the same group.
    """
    global sequence_matrix
    # Offset of 1 since we won't compare sequences to themselves
    sequence_matrix = np.full((len(sequences), len(sequences)), 0)

def read_loci():
    """
    Reads the loci file and updates the sequences appropriately.
    """
    batch = {}
    for line in open(input_file_name, 'r'):
        if line[0:2] == '//': # we've reached the end of the batch
            process_batch(batch)
            batch = {}
        else: # add the sequence to the batch
            id = line.split(' ')[0]
            sequence = line.split(' ')[-1].rstrip()
            batch[id] = sequence

def process_batch(batch):
    """
    Coordinates sequence updates based on whether or not the given batch matches the
    requirements provided in the arguments file

    Args:
        batch: A collection of id-sequence pairs read in from the loci file.
    """
    if batch_is_valid(list(batch.keys())):
        update_sequences(batch)
        update_sequence_matrix(list(batch.keys()))
    # otherwise we ignore the batch

def batch_is_valid(batch):
    """
    Determines whether or not a given group of species matches one of the desired groups.

    Args:
        batch: A list of sequence ids corresponding to a batch in the loci file.
    Returns:
        A boolean representing whether or not the batch is valid.
    """
    for group in loci_groups:
        if len(set_intersection(group['species'], batch)) >= group['minimum_matches']:
            return True
    return False

def set_intersection(list1, list2):
    """
    Performs set intersection on two given lists.
    """
    return [value for value in list1 if value in list2]

def update_sequences(batch):
    """
    Updates the sequence hash appropriately given a batch of sequences.

    Args:
        batch: A collection of id-sequence pairs.
    """
    sequence_length = len(list(batch.values())[0])
    for key, value in batch.items():
        sequences[key] += value.replace('-', 'N') # update species in the batch with their sequence
    for species in sequences.keys():
        if species not in batch.keys(): # update species not in the batch with N's
            sequences[species] += 'N' * sequence_length

def update_sequence_matrix(batch_ids):
    """
    Updates the sequence matrix given the ids of sequences in a batch.

    Args:
        batch_ids: A list of species ids
    """
    id_positions = [list(sequences.keys()).index(sequence) for sequence in batch_ids]
    # All permutations of length 2
    all_permutations = itertools.permutations(id_positions, 2)
    for position in all_permutations:
        sequence_matrix[position[0], position[1]] += 1

def write_loci_to_phylip():
    """
    Writes the sequences to the .phy output file.
    """
    output_file = open(output_file_name, 'w')
    output_file.write(str(len(sequences.keys())) + " " + str(len(list(sequences.values())[0])))
    for id, sequence in sequences.items():
        output_file.write("\n" + id + "\n" + sequence)
    output_file.close()

def write_sequence_matrix():
    """
    Writes the sequence matrix to the output file.
    """
    with open(matrix_file_name, "w") as f:
        # Write top line of species ids
        f.write('\t')
        for id in list(sequences.keys()):
            f.write(id + '\t')
        f.write('\n')

        # Write rows like <sequence id> <counts>
        index = 0
        for row in sequence_matrix:
            f.write(list(sequences.keys())[index] + '\t')
            f.write(np.array2string(row,
                                    separator='\t',
                                    max_line_width=np.nan,
                                    formatter={'int':lambda x: "%4d" % x})[1:-1] + '\n')
            index += 1

def main():
    """
    The driver function, which is called when the script is executed.
    """
    initialize_variables(sys.argv)
    read_loci()
    write_loci_to_phylip()
    write_sequence_matrix()

if __name__ == "__main__": main()
