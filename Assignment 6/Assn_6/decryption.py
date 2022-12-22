# -*- coding: utf-8 -*-

import numpy as np
import random
import sympy as sp
from numpy.linalg import matrix_rank
np.set_printoptions(threshold=np.inf)
from pyfinite import ffield


str_to_hex = {'f' : [0,0,0,0],
 'g' : [0,0,0,1],
 'h' : [0,0,1,0],
 'i' : [0,0,1,1],
 'j' : [0,1,0,0],
 'k' : [0,1,0,1],
 'l' : [0,1,1,0],
 'm' : [0,1,1,1],
 'n' : [1,0,0,0],
 'o' : [1,0,0,1],
 'p' : [1,0,1,0],
 'q' : [1,0,1,1],
 'r' : [1,1,0,0],
 's' : [1,1,0,1],
 't' : [1,1,1,0],
 'u' : [1,1,1,1]}

hex_to_str = {'0000': 'f',
 '0001': 'g',
 '0010': 'h',
 '0011': 'i',
 '0100': 'j',
 '0101': 'k',
 '0110': 'l',
 '0111': 'm',
 '1000': 'n',
 '1001': 'o',
 '1010': 'p',
 '1011': 'q',
 '1100': 'r',
 '1101': 's',
 '1110': 't',
 '1111': 'u'}

#It converts byte to corresponding two charachters
def byte_str(b):
    binnum = '{:0>8}'.format(format(b,"b"))
    a = hex_to_str[binnum[0:4]], hex_to_str[binnum[4:8]]
    return a[0]+a[1]
#It maps set of two charachters(byte) to it's corresponding hex value
def map_to_Ascii(str):
    char = chr(16*(ord(str[0]) - ord('f')) + ord(str[1]) - ord('f'))
    return char
#It takes full block of ip eg 'fffffffffffffffg' and outputs corresponding hex list
def block_to_Ascii(ch):
    plaintext = ""
    for i in range(0, len(ch), 2):
        plaintext += map_to_Ascii(ch[i:i+2])
    return plaintext

#It contains all the required functions for reuse
#Add, Multiply,Expo, add_vectors, Multiplyscalars, Lin_Trans
exp_store = [[-1]*128 for i in range(128)]

F = ffield.FField(7)

def Add (num1, num2):
    return int(num1) ^ int(num2)

def Multiply (num1, num2):
    return F.Multiply(num1, num2)

def Expo (num, pow):
    if exp_store[num][pow] != -1:
        return exp_store[num][pow]

    result = 0;
    if pow == 0:
        result = 1
    elif pow == 1:
        result = num
    elif pow%2 == 0:
        sqrt_no = Expo(num, pow>>1)
        result = Multiply(sqrt_no, sqrt_no)
    else:
        sqrt_no = Expo(num, pow>>1)
        result = Multiply(sqrt_no, sqrt_no)
        result = Multiply(num, result)

    exp_store[num][pow] = result
    return result

def add_vectors (v1, v2):
    result = [0]*8
    for i, (e1, e2) in enumerate(zip(v1, v2)):
        result[i] = Add(e1, e2)
    return result

def Multiplyscalars (v, elem):
    result = [0]*8
    for i, e in enumerate(v):
        result[i] = Multiply(e, elem)
    return result

def Lin_Trans (mat, elist):
    result = [0]*8
    for row, elem in zip(mat, elist):
        result = add_vectors(Multiplyscalars(row, elem), result)
    return result

#This list will consist of all possible exponents
exp_poss_val = [[] for i in range(8)]
#This list will consist of all possible diagonal values
diag_poss_val = [[[] for i in range(8)] for j in range(8)]
input_file = open("inputs.txt", 'r')
output_file = open("outputs.txt", 'r')
for ind, (iline, oline) in enumerate(zip(input_file.readlines(), output_file.readlines())):
    inpString = []
    outString = []
    #Converting each input to corresponding hex values
    for hexi in iline.strip().split(" "):
        inpString.append(block_to_Ascii(hexi)[ind])
    for hexo in oline.strip().split(" "):
        outString.append(block_to_Ascii(hexo)[ind])
        
    for i in range(1, 127):
        for j in range(1, 128):
            flag = True
            for inp, out in zip(inpString, outString):
                #We iterate over all possible values of ei and ajj to check if input maps to output
                if ord(out) != Expo(Multiply(Expo(Multiply(Expo(ord(inp), i), j), i), j), i):
                    flag = False
                    break
            if flag:
                #If yes, then we append them to corresponding possible lists
                exp_poss_val[ind].append(i)
                diag_poss_val[ind][ind].append(j)
print(diag_poss_val)
print(exp_poss_val)

input_file = open("inputs.txt", 'r')
output_file = open("outputs.txt", 'r')
for ind, (iline, oline) in enumerate(zip(input_file.readlines(), output_file.readlines())):
    #Considering only first 6 pairs
    if ind > 6 :
        break
    inpString = []
    outString = []
    #Converting each input to corresponding hex values
    for hexi in iline.strip().split(" "):
        inpString.append(block_to_Ascii(hexi)[ind]) 
    for hexo in oline.strip().split(" "):
        outString.append(block_to_Ascii(hexo)[ind+1])
    for i in range(1, 128):
        #We iterate over all possible pairs of exponents and diagonal values
        for p1, e1 in zip(exp_poss_val[ind+1], diag_poss_val[ind+1][ind+1]):
            for p2, e2 in zip(exp_poss_val[ind], diag_poss_val[ind][ind]):
                flag = True
                for input, output in zip(inpString, outString):
                    #We substitute the pairs ad=nd iterate over all possible values of i
                    #We do this by forming virtual triangle (aii,aij,ajj)
                    if ord(output) != Expo(Add(Multiply(Expo(Multiply(Expo(ord(input), p2), e2), p2), i) ,Multiply(Expo(Multiply(Expo(ord(input), p2), i), p1), e1)), p1):
                        flag = False
                        break
                if flag:
                    #If we find such value, then we can discard other possibilities and finalize the lists
                    exp_poss_val[ind+1] = [p1]
                    diag_poss_val[ind+1][ind+1] = [e1]
                    exp_poss_val[ind] = [p2]
                    diag_poss_val[ind][ind] = [e2]
                    diag_poss_val[ind][ind+1] = [i]
print(diag_poss_val)
print(exp_poss_val)

#Function used to apply operations (eaeae)
def eaeae (plaintext, lin_mat, exp_mat):
    plaintext = [ord(ch) for ch in plaintext]
    output = [[0 for j in range (8)] for i in range(8)]
    for ind, elem in enumerate(plaintext):
        output[0][ind] = Expo(elem, exp_mat[ind])

    output[1] = Lin_Trans(lin_mat, output[0])

    for ind, elem in enumerate(output[1]):
        output[2][ind] = Expo(elem, exp_mat[ind])

    output[3] = Lin_Trans(lin_mat, output[2])

    for ind, elem in enumerate(output[3]):
        output[4][ind] = Expo(elem, exp_mat[ind])
    return output[4]
for index in range(6):
    #As we have already found element next to diagonal thus skipping two elements every time
    of = index + 2
    
    Exp_mat = [e[0] for e in exp_poss_val]
    lin_trans_mat = [[0 for i in range(8)] for j in range(8)]
    #We fill all the empty [] elements with 0
    for i in range(8):
        for j in range(8):
            lin_trans_mat[i][j] = 0 if len(diag_poss_val[i][j]) == 0 else diag_poss_val[i][j][0]
    input = open("inputs.txt", 'r')
    out = open("outputs.txt", 'r')
    for ind, (iline, oline) in enumerate(zip(input.readlines(), out.readlines())):
        if ind > (7-of):
            continue
        inpString = [block_to_Ascii(msg) for msg in iline.strip().split(" ")]
        outString = [block_to_Ascii(msg) for msg in oline.strip().split(" ")]
        #We iterate over all possible values of ai,j to find which one satisfies eaeae = output
        for i in range(1, 128):
            lin_trans_mat[ind][ind+of] = i
            flag = True
            for inp, out in zip(inpString, outString):
                if eaeae(inp, lin_trans_mat, Exp_mat)[ind+of] != ord(out[ind+of]):
                    flag = False
                    break
            if flag:
                diag_poss_val[ind][ind+of] = [i]
    input.close();
    out.close();


lin_trans_mat = [[0 for i in range(8)] for j in range(8)]
for i in range(8):
    for j in range(8):
        lin_trans_mat[i][j] = 0 if len(diag_poss_val[i][j]) == 0 else diag_poss_val[i][j][0]

print(lin_trans_mat)
print(Exp_mat)

#Two halves of password
password_1 = "lhlrmfjffiitflmi"
password_2 = "jtlimshtmhfrfqkk"

#We iterate over all possible 128 byte values and perform eaeae to check for which input the password (half) matches
def decrypt(password):
    password_ = block_to_Ascii(password)
    res = ""
    for ind in range(8):
        for ans in range(128):
            input = res + byte_str(ans)+(16-len(res)-2)*'f'
            if ord(password_[ind]) == eaeae(block_to_Ascii(input),lin_trans_mat, Exp_mat)[ind]:
                res += byte_str(ans)
                break
    return res
decrypt_pass = block_to_Ascii(decrypt(password_1))+block_to_Ascii(decrypt(password_2)) 
print(block_to_Ascii(decrypt(password_1))+block_to_Ascii(decrypt(password_2)))

final_pass = decrypt_pass[:-6]
print("Final Passowrd is : ",final_pass)     

