import pandas as pd

def add_line(width, noMagnets, effectiveMass, muSc, alpha, M, addedSinu,ratio):
    df = pd.read_csv("data/wiresdata.csv", sep=',')

    newline = df.append({'W' : width,
                        'Ns[i]' : noMagnets,
                        'eM' : effectiveMass,
                        'mu' : muSc,
                        'al' : alpha,
                        'M' : M,
                        'int(added)' : addedSinu,
                        'ratio' : ratio}, ignore_index=True)

    print(newline)

    newline.to_csv('data/wiresdata.csv',sep=',',index=False)

    print(df)

if __name__ == "__main__":
    add_line(1,2,3,4,5,6,True,123)