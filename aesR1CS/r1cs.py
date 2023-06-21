from aes128_r1cs_const import *


def get_r1cs(constraints_count, constaint_len, witness_index_base, aes_size):
    result = {
        'A': [[0 for _ in range(constaint_len)] for _ in range(constraints_count)],
        'B': [[0 for _ in range(constaint_len)] for _ in range(constraints_count)],
        'C': [[0 for _ in range(constaint_len)] for _ in range(constraints_count)]}

    # AES round operations
    aes_output_index_base = 1
    aes_input_index_base = 17
    key_index_base = witness_index_base + aes_params[aes_size]['Nr'] * 256
    if aes_size in [192, 256]:
        aes_input_index_base = 33
        key_index_base = witness_index_base + 2 * aes_params[aes_size]['Nr'] * 256
    key_hy_index_base = key_index_base + 4 * (aes_params[aes_size]['Nk'] + aes_params[aes_size]['Nkr'])
    for i in range(aes_params[aes_size]['Nr']):
        for j in range(16):
            if i == 0:
                result['A'][j][aes_input_index_base + j] = 1
                result['A'][j][key_index_base + j] = 1
            else:
                result['A'][i * 288 + j][0] = 0x63
                mc_base = witness_index_base + 256 * (i - 1) + 128
                result['A'][i * 288 + j][mc_base:mc_base + 128] = mc[j]
                result['A'][i * 288 + j][key_index_base:key_index_base +
                                         aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][i * 16 + j]
            c_index_base = witness_index_base + i * 256 + 8 * j
            b_index_base = c_index_base + 128
            result['B'][i * 288 + j][b_index_base:b_index_base + 8] = v1
            result['C'][i * 288 + j][0] = 1
            result['C'][i * 288 + j][c_index_base:c_index_base + 8] = vf
        for j in range(256):
            if j < 128:
                inner_index_base = witness_index_base + 128
            else:
                inner_index_base = witness_index_base - 128
            result['A'][i * 288 + 16 +
                        j][inner_index_base + 256*i + j] = 1
            result['B'][i * 288 + 16 + j][0] = 1
            result['B'][i * 288 + 16 +
                        j][inner_index_base + 256*i + j] = 1

        for j in range(16):
            cur_h_index_base = witness_index_base + 256 * i + 8 * j
            cur_g_index_base = cur_h_index_base + 128
            result['A'][i * 288 + 272 +
                        j][cur_h_index_base:cur_h_index_base + 8] = e7
            result['B'][i * 288 + 272 + j][0] = 128
            result['B'][i * 288 + 272 +
                        j][cur_h_index_base:cur_h_index_base + 8] = v1
            result['B'][i * 288 + 272 +
                        j][cur_g_index_base:cur_g_index_base + 8] = v2

        if i == aes_params[aes_size]['Nr'] - 1:
            for j in range(16):
                result['A'][(i + 1) * 288 + j][0] = 0x63
                a_index_base = witness_index_base + 256*i + 128
                result['A'][(i + 1) * 288 +
                            j][a_index_base:a_index_base + 128] = ark[j]
                result['A'][(i + 1) * 288 + j][key_index_base:key_index_base +
                                               aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][(i + 1) * 16 + j]
                result['B'][(i + 1) * 288 + j][0] = 1
                result['C'][(i + 1) * 288 + j][aes_output_index_base + j] = 1

    if aes_size in [192, 256]:
        constraint_index_base = aes_params[aes_size]['Nr'] * 288 + 16
        aes_output_index_base = 17
        aes_input_index_base = 49
        # Doing 2 AES
        first_aes_data_index = witness_index_base + aes_params[aes_size]['Nr'] * 256
        for i in range(aes_params[aes_size]['Nr']):
            for j in range(16):
                if i == 0:
                    result['A'][constraint_index_base +
                                j][aes_input_index_base + j] = 1
                    result['A'][constraint_index_base +
                                j][key_index_base + j] = 1
                else:
                    result['A'][constraint_index_base + i * 288 + j][0] = 0x63
                    mc_base = first_aes_data_index + 256 * (i - 1) + 128
                    result['A'][constraint_index_base + i *
                                288 + j][mc_base:mc_base + 128] = mc[j]
                    result['A'][constraint_index_base + i * 288 + j][key_index_base:key_index_base +
                                                                     aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][i * 16 + j]
                c_index_base = first_aes_data_index + i * 256 + 8 * j
                b_index_base = c_index_base + 128
                result['B'][constraint_index_base + i * 288 +
                            j][b_index_base:b_index_base + 8] = v1
                result['C'][constraint_index_base + i * 288 + j][0] = 1
                result['C'][constraint_index_base + i * 288 +
                            j][c_index_base:c_index_base + 8] = vf
            for j in range(256):
                if j < 128:
                    inner_index_base = first_aes_data_index + 128
                else:
                    inner_index_base = first_aes_data_index - 128
                result['A'][constraint_index_base + i * 288 + 16 +
                            j][inner_index_base + 256*i + j] = 1
                result['B'][constraint_index_base + i * 288 + 16 + j][0] = 1
                result['B'][constraint_index_base + i * 288 + 16 +
                            j][inner_index_base + 256*i + j] = 1

            for j in range(16):
                cur_h_index_base = first_aes_data_index + 256 * i + 8 * j
                cur_g_index_base = cur_h_index_base + 128
                result['A'][constraint_index_base + i * 288 + 272 +
                            j][cur_h_index_base:cur_h_index_base + 8] = e7
                result['B'][constraint_index_base + i * 288 + 272 + j][0] = 128
                result['B'][constraint_index_base + i * 288 + 272 +
                            j][cur_h_index_base:cur_h_index_base + 8] = v1
                result['B'][constraint_index_base + i * 288 + 272 +
                            j][cur_g_index_base:cur_g_index_base + 8] = v2

            if i == aes_params[aes_size]['Nr'] - 1:
                for j in range(16):
                    result['A'][constraint_index_base +
                                (i + 1) * 288 + j][0] = 0x63
                    a_index_base = first_aes_data_index + 256*i + 128
                    result['A'][constraint_index_base + (i + 1) * 288 +
                                j][a_index_base:a_index_base + 128] = ark[j]
                    ak_index_base = first_aes_data_index + \
                        aes_params[aes_size]['Nr'] * 256
                    result['A'][constraint_index_base + (i + 1) * 288 + j][ak_index_base:ak_index_base +
                                                                           aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][(i + 1) * 16 + j]
                    result['B'][constraint_index_base +
                                (i + 1) * 288 + j][0] = 1
                    result['C'][constraint_index_base + (i + 1) * 288 +
                                j][aes_output_index_base + j] = 1

    # 288 * aes_params[aes_size]['Nr'] + 16 constraints added
    # AES key schedule
    ks_constraint_index = aes_params[aes_size]['Nr'] * 288 + 16
    if aes_size in [192, 256]:
        ks_constraint_index = (aes_params[aes_size]['Nr'] * 288 + 16) * 2
    for i in range(aes_params[aes_size]['Nkr']):

        for j in range(4):
            # Iterating R_{i,j}
            result['A'][ks_constraint_index + 76 * i + j][key_index_base:key_index_base +
                                                          aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][aes_params[aes_size]['key_schedule_index'](i, j)]
            c_index_base = key_hy_index_base + i*64 + 8 * ((j + 1) % 4)
            if aes_size == 256 and i % 2:
                c_index_base = key_hy_index_base + i*64 + 8 * j
            b_index_base = c_index_base + 32
            result['B'][ks_constraint_index + 76 * i +
                        j][b_index_base:b_index_base + 8] = v1
            result['C'][ks_constraint_index + 76 * i + j][0] = 1
            result['C'][ks_constraint_index + 76 * i +
                        j][c_index_base:c_index_base + 8] = vf
        for j in range(64):
            if j < 32:
                inner_index_base = key_hy_index_base + 32
            else:
                inner_index_base = key_hy_index_base - 32
            result['A'][ks_constraint_index + 76 * i + 4 +
                        j][inner_index_base + i*64 + j] = 1
            result['B'][ks_constraint_index + 76 * i + 4 + j][0] = 1
            result['B'][ks_constraint_index + 76 * i +
                        4 + j][inner_index_base + i*64 + j] = 1

        for j in range(4):
            cur_kh_index_base = key_hy_index_base + i * 64 + 8 * j
            cur_kg_index_base = cur_kh_index_base + 32
            result['A'][ks_constraint_index + 76 * i + 68 +
                        j][cur_kh_index_base:cur_kh_index_base + 8] = e7
            result['B'][ks_constraint_index + 76 * i + 68 + j][0] = 128
            result['B'][ks_constraint_index + 76 * i + 68 +
                        j][cur_kh_index_base:cur_kh_index_base + 8] = v1
            result['B'][ks_constraint_index + 76 * i + 68 +
                        j][cur_kg_index_base:cur_kg_index_base + 8] = v2
        for j in range(4):
            if j == 0:
                if aes_size == 256:
                    if i % 2:
                        result['A'][ks_constraint_index + 76 *
                                    i + 72 + j][0] = 0x63
                    else:
                        result['A'][ks_constraint_index + 76 *
                                    i + 72 + j][0] = 0x63 ^ rci[i // 2]
                else:
                    result['A'][ks_constraint_index + 76 *
                                i + 72 + j][0] = 0x63 ^ rci[i]
            else:
                result['A'][ks_constraint_index + 76 * i + 72 + j][0] = 0x63
            a_index_base = key_hy_index_base + i*64 + 32 + 8 * ((j + 1) % 4)
            if aes_size == 256 and i % 2:
                a_index_base = key_hy_index_base + i*64 + 32 + 8 * j
            result['A'][ks_constraint_index + 76 * i +
                        72 + j][a_index_base:a_index_base + 8] = w1
            result['B'][ks_constraint_index + 76 * i + 72 + j][0] = 1
            result['C'][ks_constraint_index + 76 * i + 72 +
                        j][key_index_base + 4 * aes_params[aes_size]['Nk'] + 4 * i + j] = 1

    return result
















if __name__ == '__main__' :
    print( get_r1cs(4096, 4096, 64, 128) )