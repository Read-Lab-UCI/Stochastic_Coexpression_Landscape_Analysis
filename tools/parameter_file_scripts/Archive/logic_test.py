import itertools
g_names_a = ['g0', 'g1a']
g_names_b = ['g0', 'g1b']
header_1 = ['', 'A_00', 'A_01', 'A_10', 'A_11', '', 'B_00', 'B_01', 'B_10', 'B_11']
header_2 = ['', 'BA', 'BA', 'BA', 'BA', '', 'AB', 'AB', 'AB', 'AB']

variedLogicMatrix_a = list(itertools.product(g_names_a, repeat=4))
variedLogicMatrix_b = list(itertools.product(g_names_b, repeat=4))

print('  | '.join(str(x).rjust(4) for x in header_2))
print('  | '.join(str(x).rjust(4) for x in header_1))

for i, row_a in enumerate(variedLogicMatrix_a):
    row_b = variedLogicMatrix_b[i]
    mod_row = (i+1,) + row_a + ('',) + row_b
    print('  | '.join(str(x).rjust(4) for x in mod_row))
