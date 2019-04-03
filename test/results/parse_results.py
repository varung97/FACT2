import sys
import csv
import numpy as np
import matplotlib.pyplot as plt


def parse_and_draw(xlabel, scen_num):
    input_filename = 'data/varying' + xlabel + scen_num + '.txt'
    graph_path = 'images/varying' + xlabel + scen_num
    latex_path = 'tex/varying' + xlabel + scen_num

    with open(input_filename) as inf:
        # data ordering: ['new_w', 'new_fc', 'old_w_kn2', 'old_w_k2n', 'old_fc']

        data = inf.read().splitlines()
        data_table = []
        num_trials = int(data[0])
        data_row = 1
        while data_row < len(data):
            n, k = data[data_row].split(', ')
            data_row += 1

            times = np.array(data[data_row : data_row + 5 * num_trials], dtype=float).reshape(num_trials, 5)

            data_table.append(np.c_[np.full(num_trials, n, dtype=float), np.full(num_trials, k, dtype=float), times])

            data_row = data_row + 5 * num_trials

        data_table = np.array(data_table).swapaxes(0, 1)
        averaged_data = np.average(data_table, axis=0)

        if xlabel == 'n':
            x = averaged_data[:, 0]
        else:
            x = averaged_data[:, 1]

        plt.figure(figsize=(7, 6))
        new_w, old_w_kn2, old_w_k2n = plt.plot(x, averaged_data[:, 2], 'r--', x, averaged_data[:, 4], 'b--', x, averaged_data[:, 5], 'g--')
        plt.xlabel(xlabel)
        plt.ylabel('Time / s')
        plt.legend((new_w, old_w_kn2, old_w_k2n), ('$kn\,log\,n$', '$kn^2$', '$k^2n$'), loc='upper left')
        plt.savefig(graph_path + '_weighting.png', dpi=300, bbox_inches='tight')

        with open(latex_path + '_weighting.txt', 'w') as f:
            latex_data = [[str(int(x[i]))] + ['{:.2f}'.format(j) for j in row] for i, row in enumerate(averaged_data[:, (2, 4, 5)])]
            f.write('\\\\\n'.join([' & '.join(row) for row in latex_data]))
            f.write('\\\\\n')

        plt.figure(figsize=(7, 6))
        new_fc, old_fc = plt.plot(x, averaged_data[:, 3], 'r--', x, averaged_data[:, 6], 'b--')
        plt.xlabel(xlabel)
        plt.ylabel('Time / s')
        plt.legend((new_fc, old_fc), ('$n\,log\,n$', '$n\,log^2n$'), loc='upper left')
        plt.savefig(graph_path + '_filter.png', dpi=300, bbox_inches='tight')

        with open(latex_path + '_filter.txt', 'w') as f:
            latex_data = [[str(int(x[i]))] + ['{:.2f}'.format(j) for j in row] for i, row in enumerate(averaged_data[:, (3, 6)])]
            f.write('\\\\\n'.join([' & '.join(row) for row in latex_data]))
            f.write('\\\\\n')

        # print(np.std(data_table, axis=0) / np.average(data_table, axis=0))



if __name__ == '__main__':
    np.set_printoptions(formatter={'float': '{: 0.2f}'.format}, suppress=True)

    if sys.argv[1] == 'all':
        parse_and_draw('n', '1')
        parse_and_draw('n', '2')
        parse_and_draw('k', '1')
        parse_and_draw('k', '2')
    else:
        parse_and_draw(xlabel=sys.argv[1], scen_num=sys.argv[2])
