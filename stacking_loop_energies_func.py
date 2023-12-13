from free_energy import stacking_energies
# moved to free_energy.py
def eS(i, j, RNA_seq, backtrace):
    size = j - i
    # find out if (i,j) closes a stacking loop?
    ext_closing = str(RNA_seq[i]) + str(RNA_seq[j])
    ext_closing = ''.join(sorted(ext_closing))

    int_clos_i = i + 1
    int_clos_j = j - 1

    int_closing = str(RNA_seq[int_clos_i]) + str(RNA_seq[int_clos_j])
    int_closing = ''.join(sorted(int_closing))

    energy = stacking_energies[ext_closing][int_closing]
    return energy/10.0    

