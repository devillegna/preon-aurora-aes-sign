from aes128_r1cs_const import *

class sp_mat:
    def __init__(self,n_cols,n_rows):
        self.n_cols = n_cols
        self.n_rows = n_rows
        self.col = [ [] for i in range(n_cols) ]

def _remove_0( list_of_i_v ):
    r = []
    for i_v in list_of_i_v :
        if 0 != i_v[1] : r.append(i_v)
    return r

def get_aes128r1cs(constraints_count, constaint_len, witness_index_base, aes_size):
    mat_a = sp_mat(constaint_len,constraints_count)
    mat_b = sp_mat(constaint_len,constraints_count)
    mat_c = sp_mat(constaint_len,constraints_count)

    #result = {
    #    'A': [[0 for _ in range(constaint_len)] for _ in range(constraints_count)],
    #    'B': [[0 for _ in range(constaint_len)] for _ in range(constraints_count)],
    #    'C': [[0 for _ in range(constaint_len)] for _ in range(constraints_count)]}

    # AES round operations
    aes_output_index_base = 1
    aes_input_index_base = 17
    key_index_base = witness_index_base + aes_params[aes_size]['Nr'] * 256
    key_hy_index_base = key_index_base + 4 * (aes_params[aes_size]['Nk'] + aes_params[aes_size]['Nkr'])

    for i in range(aes_params[aes_size]['Nr']):
        for j in range(16):
            if i == 0:
                #result['A'][j][aes_input_index_base + j] = 1
                #result['A'][j][key_index_base + j] = 1
                mat_a.col[aes_input_index_base + j].append( (j,1) )
                mat_a.col[key_index_base       + j].append( (j,1) )
            else:
                #result['A'][i * 288 + j][0] = 0x63
                mc_base = witness_index_base + 256 * (i - 1) + 128
                #result['A'][i * 288 + j][mc_base:mc_base + 128] = mc[j]
                #result['A'][i * 288 + j][key_index_base:key_index_base +
                #                         aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][i * 16 + j]
                mat_a.col[0].append( (i*288+j,0x63) )
                for k in range(128): mat_a.col[mc_base+k].append( (i*288+j,mc[j][k]) )
                for k in range(aes_params[aes_size]['key_basis_len']) :
                    mat_a.col[key_index_base+k].append( (i*288+j,aes_params[aes_size]['key_basis'][i*16+j][k]) )
            c_index_base = witness_index_base + i * 256 + 8 * j
            b_index_base = c_index_base + 128
            #result['B'][i * 288 + j][b_index_base:b_index_base + 8] = v1
            #result['C'][i * 288 + j][0] = 1
            #result['C'][i * 288 + j][c_index_base:c_index_base + 8] = vf
            for k in range(8): mat_b.col[b_index_base+k].append( (i*288+j,v1[k]) )
            mat_c.col[0].append( (i*288+j,1) )
            for k in range(8): mat_c.col[c_index_base+k].append( (i*288+j,vf[k]) )
        for j in range(256):
            if j < 128:
                inner_index_base = witness_index_base + 128
            else:
                inner_index_base = witness_index_base - 128
            #result['A'][i * 288 + 16 + j][inner_index_base + 256*i + j] = 1
            #result['B'][i * 288 + 16 + j][0] = 1
            #result['B'][i * 288 + 16 + j][inner_index_base + 256*i + j] = 1
            mat_a.col[inner_index_base + 256*i + j].append( (i*288+16+j,1) )
            mat_b.col[                           0].append( (i*288+16+j,1) )
            mat_b.col[inner_index_base + 256*i + j].append( (i*288+16+j,1) )
        for j in range(16):
            cur_h_index_base = witness_index_base + 256 * i + 8 * j
            cur_g_index_base = cur_h_index_base + 128
            #result['A'][i * 288 + 272 + j][cur_h_index_base:cur_h_index_base + 8] = e7
            #result['B'][i * 288 + 272 + j][0] = 128
            #result['B'][i * 288 + 272 + j][cur_h_index_base:cur_h_index_base + 8] = v1
            #result['B'][i * 288 + 272 + j][cur_g_index_base:cur_g_index_base + 8] = v2
            for k in range(8): mat_a.col[cur_h_index_base+k].append( (i*288+272+j,e7[k]) )
            mat_b.col[0].append( (i*288+272+j,128) )
            for k in range(8): mat_b.col[cur_h_index_base+k].append( (i*288+272+j,v1[k]) )
            for k in range(8): mat_b.col[cur_g_index_base+k].append( (i*288+272+j,v2[k]) )

        if i == aes_params[aes_size]['Nr'] - 1:
            for j in range(16):
                #result['A'][(i + 1) * 288 + j][0] = 0x63
                a_index_base = witness_index_base + 256*i + 128
                #result['A'][(i + 1) * 288 + j][a_index_base:a_index_base + 128] = ark[j]
                #result['A'][(i + 1) * 288 + j][key_index_base:key_index_base + aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][(i + 1) * 16 + j]
                #result['B'][(i + 1) * 288 + j][0] = 1
                #result['C'][(i + 1) * 288 + j][aes_output_index_base + j] = 1
                mat_a.col[0].append( (i*288+288+j,0x63) )
                for k in range(128): mat_a.col[a_index_base+k].append( (i*288+288+j,ark[j][k]) )
                for k in range(aes_params[aes_size]['key_basis_len']):
                    mat_a.col[key_index_base+k].append( (i*288+288+j,aes_params[aes_size]['key_basis'][(i+1)*16+j][k]) )
                mat_b.col[                        0].append( (i*288+288+j,1) )
                mat_c.col[aes_output_index_base + j].append( (i*288+288+j,1) )
    # 288 * aes_params[aes_size]['Nr'] + 16 constraints added
    # AES key schedule
    ks_constraint_index = aes_params[aes_size]['Nr'] * 288 + 16
    for i in range(aes_params[aes_size]['Nkr']):
        for j in range(4):
            # Iterating R_{i,j}
            #result['A'][ks_constraint_index + 76 * i + j][key_index_base:key_index_base +
            #                                              aes_params[aes_size]['key_basis_len']] = aes_params[aes_size]['key_basis'][aes_params[aes_size]['key_schedule_index'](i, j)]
            for k in range(aes_params[aes_size]['key_basis_len']): 
                mat_a.col[key_index_base+k].append( 
                (ks_constraint_index+76*i+j,aes_params[aes_size]['key_basis'][aes_params[aes_size]['key_schedule_index'](i, j)][k]) )
            c_index_base = key_hy_index_base + i*64 + 8 * ((j + 1) % 4)
            b_index_base = c_index_base + 32
            #result['B'][ks_constraint_index + 76 * i + j][b_index_base:b_index_base + 8] = v1
            #result['C'][ks_constraint_index + 76 * i + j][0] = 1
            #result['C'][ks_constraint_index + 76 * i + j][c_index_base:c_index_base + 8] = vf
            for k in range(8): mat_b.col[b_index_base+k].append( (ks_constraint_index+76*i+j,v1[k]) )
            mat_c.col[0].append( (ks_constraint_index + 76 * i + j,1) )
            for k in range(8): mat_c.col[c_index_base+k].append( (ks_constraint_index + 76 * i + j,vf[k]) )
        for j in range(64):
            if j < 32:
                inner_index_base = key_hy_index_base + 32
            else:
                inner_index_base = key_hy_index_base - 32
            #result['A'][ks_constraint_index + 76 * i + 4 + j][inner_index_base + i*64 + j] = 1
            #result['B'][ks_constraint_index + 76 * i + 4 + j][0] = 1
            #result['B'][ks_constraint_index + 76 * i + 4 + j][inner_index_base + i*64 + j] = 1
            mat_a.col[inner_index_base + i*64 + j].append( (ks_constraint_index + 76 * i + 4 + j,1) )
            mat_b.col[                          0].append( (ks_constraint_index + 76 * i + 4 + j,1) )
            mat_b.col[inner_index_base + i*64 + j].append( (ks_constraint_index + 76 * i + 4 + j,1) )
        for j in range(4):
            cur_kh_index_base = key_hy_index_base + i * 64 + 8 * j
            cur_kg_index_base = cur_kh_index_base + 32
            #result['A'][ks_constraint_index + 76 * i + 68 + j][cur_kh_index_base:cur_kh_index_base + 8] = e7
            #result['B'][ks_constraint_index + 76 * i + 68 + j][0] = 128
            #result['B'][ks_constraint_index + 76 * i + 68 + j][cur_kh_index_base:cur_kh_index_base + 8] = v1
            #result['B'][ks_constraint_index + 76 * i + 68 + j][cur_kg_index_base:cur_kg_index_base + 8] = v2
            for k in range(8): mat_a.col[cur_kh_index_base+k].append( (ks_constraint_index+76*i+68+j,e7[k]) )
            mat_b.col[                                     0].append( (ks_constraint_index+76*i+68+j,128) )
            for k in range(8): mat_b.col[cur_kh_index_base+k].append( (ks_constraint_index+76*i+68+j,v1[k]) )
            for k in range(8): mat_b.col[cur_kg_index_base+k].append( (ks_constraint_index+76*i+68+j,v2[k]) )
        for j in range(4):
            if j == 0:
                #result['A'][ks_constraint_index + 76 * i + 72 + j][0] = 0x63 ^ rci[i]
                mat_a.col[0].append( (ks_constraint_index + 76 * i + 72 + j,0x63^rci[i]) )
            else:
                #result['A'][ks_constraint_index + 76 * i + 72 + j][0] = 0x63
                mat_a.col[0].append( (ks_constraint_index + 76 * i + 72 + j,0x63) )
            a_index_base = key_hy_index_base + i*64 + 32 + 8 * ((j + 1) % 4)
            #result['A'][ks_constraint_index + 76 * i + 72 + j][a_index_base:a_index_base + 8] = w1
            #result['B'][ks_constraint_index + 76 * i + 72 + j][0] = 1
            #result['C'][ks_constraint_index + 76 * i + 72 + j][key_index_base + 4 * aes_params[aes_size]['Nk'] + 4 * i + j] = 1
            for k in range(8): mat_a.col[a_index_base+k].append( (ks_constraint_index+76*i+72+j,w1[k]) )
            mat_b.col[0].append( (ks_constraint_index + 76 * i + 72 + j,1) )
            mat_c.col[key_index_base + 4 * aes_params[aes_size]['Nk'] + 4 * i + j].append( (ks_constraint_index + 76 * i + 72 + j,1) )
    #for i,c in enumerate(mat_a.col) : mat_a.col[i]=list(filter(lambda x,y:y!=0,c))
    #for i,c in enumerate(mat_b.col) : mat_b.col[i]=list(filter(lambda x,y:y!=0,c))
    #for i,c in enumerate(mat_c.col) : mat_c.col[i]=list(filter(lambda x,y:y!=0,c))
    for i,c in enumerate(mat_a.col) : mat_a.col[i]=_remove_0(c)
    for i,c in enumerate(mat_b.col) : mat_b.col[i]=_remove_0(c)
    for i,c in enumerate(mat_c.col) : mat_c.col[i]=_remove_0(c)
    return mat_a, mat_b, mat_c


