import sys
import csv


NUM_TRIALS = 3


if __name__ == '__main__':
    input_filename = sys.argv[1]
    output_filename = input_filename.strip('.txt') + '.csv'
    with open(input_filename) as inf, open(output_filename, 'w') as outf:
        csvwriter = csv.writer(outf, delimiter=',')
        csvwriter.writerow(['n', 'k', 'new_w', 'old_w_kn2', 'old_w_k2n', 'new_fc', 'old_fc'])

        data = inf.read().splitlines()
        data_row = 0
        while data_row < len(data):
            n, k = data[data_row].split(', ')
            data_row += 1
            for j in range(NUM_TRIALS):
                new_w, new_fc, old_w_kn2, old_w_k2n, old_fc = data[data_row : data_row + 5]
                csvwriter.writerow([n, k, new_w, old_w_kn2, old_w_k2n, new_fc, old_fc])
                data_row += 5
