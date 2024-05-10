import numpy as np
import sys
from Bio.Align import substitution_matrices

amino_acids = "ACDEFGHIKLMNPQRSTVWY-"

def convertdict(mat, string):
    scor_mat = {}
    for i in range(len(string)):
        for j in range(len(string)):
            scor_mat[(string[i], string[j])] = mat[string[i], string[j]]
    return scor_mat

def create_scoring_matrix(m, n, sigma):
    # Make all nodes 0 b/c local alignment can start from anywhere
    return np.zeros((m+1, n+1), dtype=int)

def create_profile(strings, symbols):
    prof = np.zeros((len(symbols), len(strings[0])))
    for string in strings:
        for i in range(len(string)):
            prof[symbols.index(string[i])][i] += 1
    #norm = np.sum(prof, axis = 0)[0]
    return prof 

def calculate_match(x, y, scoring_matrix, i, j, sigma):
    total = 0
    x = x/np.sum(x, axis = 0)[0]
    y = y/np.sum(y, axis = 0)[0]
    for k in range(len(x)):
        for l in range(len(y)):
            if amino_acids[k] == "-" or amino_acids[l] == "-":
                score = sigma
            else:
                score = scoring_matrix[(amino_acids[k], amino_acids[l])] 
            #Probably add 1 back to i and j. Quick fix for index error but we don't know why.
            total += score * x[k][i] * y[l][j]
    return total

def fill_scoring_matrix(x, y, S, sigma, scoring_matrix):
    # Fill scoring matrix
    m, n = len(x[0]), len(y[0])
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m+1):
        for j in range(1, n+1):
            #scores for possible alignments
            insert = S[i, j-1] - sigma
            delete = S[i-1, j] - sigma
            match = S[i-1, j-1] + calculate_match(x, y, scoring_matrix, i-1, j-1, sigma)

            S[i, j] = max(match, delete, insert, 0) #recurrence relation

            #keep track
            if S[i, j] > max_score:
                max_score = S[i, j]
                max_pos = (i, j)

    return S, max_score, max_pos


def backtrack(S, x, y, start_pos, sigma, scoring_matrix):
    #Backtrack to find optimal local alignment
    i, j = start_pos  # Start backtracking from the position of the max score
    i -= 1
    j -= 1
    # Backtrack until reaching a cell with score 0 (sigals end of local alignment)
    #xcopy = np.zeros((len(x), len(x[0])))
    #ycopy = np.zeros((len(y), len(y[0])))
    xcopy = np.zeros((len(amino_acids), 0))
    ycopy = np.zeros((len(amino_acids), 0))
    while i > 0 and j > 0 and S[i, j] > 0:
        if S[i, j] == S[i-1, j-1] + calculate_match(x, y, scoring_matrix, i-1, j-1, sigma):
            # Diagonal: match/mismatch
            #X_align = x[i-1] + X_align
            #Y_align = y[j-1] + Y_align
            xcopy = np.insert(xcopy, 0, x[:,i], axis=1)
            ycopy = np.insert(ycopy, 0, y[:,j], axis=1)
            i -= 1
            j -= 1
        elif S[i, j] == S[i, j-1] - sigma:
            #Left: insertion
            ins = np.zeros(len(amino_acids))
            ins[-1] = np.sum(x, axis = 0)[0]
            xcopy = np.insert(xcopy, 0, ins, axis=1)
            ycopy = np.insert(ycopy, 0, y[:,j], axis=1)
            j -= 1
        else: # S[i, j] == S[i-1, j] - sigma:
            #Up: deletion
            # X_align = x[i-1] + X_align
            # Y_align = "-" + Y_align
            ins = np.zeros(len(amino_acids))
            ins[-1] = np.sum(y, axis = 0)[0]
            xcopy = np.insert(xcopy, 0, x[:,i], axis=1)
            ycopy = np.insert(ycopy, 0, ins, axis=1)
            # xcopy[:, i] = x[:, i]
            # ycopy = np.insert(y, j - 1, ins, axis=1)
            i -= 1
    final = xcopy + ycopy
    return final


def local_alignment(x, y, sigma=5):
    #Perform local alignment using PAM250 scoring matrix
    scoring_matrix = substitution_matrices.load("PAM250")
    scoring_matrix = convertdict(scoring_matrix, amino_acids[:-1])
    m, n = len(x[0]), len(y[0])

    S = create_scoring_matrix(m, n, sigma)
    S, max_score, max_pos = fill_scoring_matrix(x, y, S, sigma, scoring_matrix)

    final = backtrack(S, x, y, max_pos, sigma, scoring_matrix)


    return max_score, final, x, y




def pairwiseAlign(genomes):
    #Aligns each genome with all of the other genomes. Puts them in list rank and then sorts
    #based on the max score. List is returned with elements of [score, xalign, yalign]
    genomeprof = []
    for seq in genomes:
        genomeprof.append(create_profile([seq], amino_acids))
    while len(genomeprof) > 1: 
        rank = []    
        for i in range((len(genomeprof))):
            for j in range(i+1, len(genomeprof)):
                align = local_alignment(genomeprof[i], genomeprof[j], 5)
                rank.append(align)
        rank = sorted(rank, key=lambda x : x[0], reverse = True)
        x, y = rank[0][2], rank[0][3]
        genomeprof2 = [arr for arr in genomeprof if not (arr.shape == x.shape and np.equal(x, arr).all() or arr.shape == y.shape and np.equal(y, arr).all())]
        #genomeprof.pop(genomeprof.index(y))
        genomeprof = genomeprof2
        genomeprof.append(rank[0][1])
    return genomeprof[0]


#Not sure if we still need this function? See create_profile()
def createProfile(alignment):
    x = alignment[1]
    y = alignment[2]
    prof = np.zeros((len(amino_acids), len(x)))
    for i in range(len(x)):
        if x[i] == y[i]:
            prof[amino_acids.index(x[i])][i] += 1
    
    return prof
    






scoring_matrix = substitution_matrices.load("PAM250")
seq = ["MKTIIALSYIFCLVCLVCLVCLVCLV","TIIALSYIFCLVCLVCLVCLVFA","ALSYIFCLVCLVCLVCLVFADYK","CLVCLVCLVCLVFADYKDDDDK","IFCLVCLVCLVFADY","SYIFCLVCLVCLVCLVFA"]
#seq = ["MKKIKIVPLILIVVVVGFGIYFYASKDKEINNTIDAIEDKNFKQVYKDSSYISKSDNGEVEMTERPIKIYNSLGVKDINIQDRKIKKVSKNKKRVDAQYKIKTNYGNIDRNVQFNFVKEDGMWKLDWDHSVIIPGMQKDQSIHIENLKSERGKILDRNNVELANTGTHMRLGIVPKNVSKKDYKAIAKELSISEDYINNKWIKIGYKMIPSFHFKTVKKMDEYLSDFAKKFHLTTNETESRNYPLEKATSHLLGYVGPINSEELKQKEYKGYKDDAVIGKKGLEKLYDKKLQHEDGYRVTIVDDNSNTIAHTLIEKKKKDGKDIQLTIDAKVQKSIYNNMKNDYGSGTAIHPQTGELLALVSTPSYDVYPFMYGMSNEEYNKLTEDKKEPLLNKFQITTSPGSTQKILTAMIGLNNKTLDDKTSYKIDGKGWQKDKSWGGYNVTRYEVVNGNIDLKQAIESSDNIFFARVALELGSKKFEKGMKKLGVGEDIPSDYPFYNAQISNKNLDNEILLADSGYGQGEILINPVQILSIYSALENNGNINAPHLLKDTKNKVWKKNIISKENINLLNDGMQQVVNKTHKEDIYRSYANLIGKSGTAELKMKQGESGRQIGWFISYDKDNPNMMMAINVKDVQDKGMASYNAKISGKVYDELYENGNKKYDIDE",
#       "MAIRIFAILFSIFSLATFAHAQEGTLERSDWRKFFSEFQAKGTIVVADERQADRAMLVFDPVRSKKRYSPASTFKIPHTLFALDAGAVRDEFQIFRWDGVNRGFAGHNQDQDLRSAMRNSTVWVYELFAKEIGDDKARRYLKKIDYGNADPSTSNGDYWIEGSLAISAQEQIAFLRKLYRNELPFRVEHQRLVKDLMIVEAGRNWILRAKTGWEGRMGWWVGWVEWPTGSVFFALNIDTPNRMDDLFKREAIVRAILRSIEALPPNPAVNSDAAR",
#       "MRNRGFGRRELLVAMAMLVSVTGCARHASGARPASTTLPAGADLADRFAELERRYDARLGVYVPATGTTAAIEYRADERFAFCSTFKAPLVAAVLHQNPLTHLDKLITYTSDDIRSISPVAQQHVQTGMTIGQLCDAAIRYSDGTAANLLLADLGGPGGGTAAFTGYLRSLGDTVSRLDAEEPELNRDPPGDERDTTTPHAIALVLQQLVLGNALPPDKRALLTDWMARNTTGAKRIRAGFPADWKVIDKTGTGDYGRANDIAVVWSPTGVPYVVAVMSDRAGGGYDAEPREALLAEAATCVAGVLA"]
# strings1 = ["ACGTCAG", "TCAGTCG"]
# strings2 = ["ATGCGCT", "CTGCGCT"]
# strings3 = [""]
# prof1 = create_profile(strings1, "ACDEFGHIKLMNPQRSTVWY-")
# prof2 = create_profile(strings2, "ACDEFGHIKLMNPQRSTVWY-")


#print(pairwiseAlign(seq))
#print({scoring_matrix})



#print(convertdict(scoring_matrix, amino_acids[:-1]))

def convertprof(filename):
    with open(filename, "r") as file:
        prof = np.loadtxt(file)
    fin = prof/np.sum(prof, axis = 0)[0]
    max = fin.max(axis = 0)
    idx = np.argmax(prof, axis = 0)
    string = ''
    seq = ''
    for i in range(len(max)):
        if max[i] == 1:
            string += "f"
            seq += amino_acids[idx[i]]
        elif max[i] >= 0.75:
            string += "v"
            seq += amino_acids[idx[i]]
        elif max[i] >= 0.5:
            string += "c"
            seq += amino_acids[idx[i]]
        else:
            string += "-"
            seq += "-"
    return string, seq

print(convertprof("bacterial_alignment_bla2.txt"))

seq_bla = ["MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDER",
           "ISLLATLPLAVHASPQPLEQIKLSESQLSGRVGMIEMDLASGRTLTAWRADERFPMMSTFKVVLCGAVLARVDAGDEQLERKIHYRQQDLVDYSPVSEKHLADGMTVGELCAAAITMSDNSAANLLLATVGGPAGLTAFLRQIGDNVTRLDRWETELNEALPGDARDTTTPASMAATLRKLLTSQRLSARSQRQLLQWMVDDRVAGPLIRSVLPAGWFIADKTGAGERGARGIVALLGPNNKAERIVVIYLRDTPASMAERN",
           "MVKKSLRQFTLMATATVTLLLGSVPLYAQTADVQQKLAELERQSGGRLGVALINTADNSQILYRADERFAMCSTSKVMAVAAVLKKSESEPNLLNQRVEIKKSDLVNYNPIAEKHVDGTMSLAELSAAALQYSDNVAMNKLISHVGGPASVTAFARQLGDETFRLDRTEPTLNTAIPGDPRDTTSPRAMAQTLRNLTLGKALGDSQRAQLVTWMKGNTTGAASIQAGLPASWVVGDKTGSGDYGTTNDIAVIWPKDRAPLILVTYFTQPQPKAESRRDVLASAAKIVTDGL",
           "FAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIA",
           "LPLAVHASPQPLEQIKLSESQLSGRVGMIEMDLASGRTLTAWRADERFPMMSTFKVVLCGAVLARVDAGDEQLERKIHYRQQDLVDYSPVSEKHLADGMTVGELCAAAITMSDNSAANLLLATVGGPAGLTAFLRQIGDNVTRLDRWETELNEALPGDARDTTTPASMAATLRKLLTSQRLSARSQRQLLQWMVDDRVAGPLIRSVLPAGWFIADKTGAGERGARGIVALLGPNNKAERIVVIYLRDTPASMAE",
           "MSIQHFRVALIPFFAAFCFPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYMTGSQATMDERNRQIAEIGASLIKHW"]



# with open("bacterial_alignment_bla2.txt", "w") as file:
#     alignment = pairwiseAlign(seq_bla)
#     for row in alignment:
#         file.write(" ".join(map(str, row)) + "\n")