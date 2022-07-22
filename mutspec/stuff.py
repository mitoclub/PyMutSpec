from mutspec.utils import CodonAnnotation

for i in range(1, 40):
    try:
        coda = CodonAnnotation(i)
        print(i, len(coda.possible_ff_contexts), len(coda.possible_syn_contexts))
    except:
        pass

# 1 48 64
# 2 48 64
# 3 48 64
# 4 48 64
# 5 48 64
# 6 48 64
# 9 48 64
# 10 48 64
# 11 48 64
# 12 48 64
# 13 48 64
# 14 48 64
# 15 48 64
# 16 48 64
# 21 48 64
# 22 48 64
# 23 48 64
# 24 48 64
# 25 48 64
# 26 48 64
# 27 48 64
# 28 48 64
# 29 64 64
# 30 48 64
# 31 48 64
# 32 48 64
# 33 48 64