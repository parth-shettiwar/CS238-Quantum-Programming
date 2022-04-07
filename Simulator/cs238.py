# -*- coding: utf-8 -*-

import numpy as np
import random
import time
import os
import sys
import itertools as itert
import cirq
import math
from pathlib import Path
import matplotlib.pyplot as plt
from cirq.contrib.qasm_import import circuit_from_qasm


curr  =-1e10

#The implementation for all gates is similar
#Implementing XGate.
def Xgate(state,m):
  n = int(np.log2(len(state))) #get the number of qubits
  new_state = np.zeros(len(state),dtype=complex) #create a new state complex vector 
  for i in range(len(state)): #iterate over all non zero entries and do the necessary operation by converting to binary form
    if(state[i]!=0):
      temp = state[i]
      state[i] = 0
      val = bin(i)[2:].zfill(n)
      bin_val = np.asarray(list(map(int,val)))
      bin_val[m] = -bin_val[m]+1
      new_val = int("".join(map(str,bin_val)),2)
      new_state[new_val] += temp
  return new_state

#Implementing Hadamard Gate  
def Hgate(state,m):
  n = int(np.log2(len(state)))
  new_state = np.zeros(len(state),dtype=complex)
  for i in range(len(state)):
    if(state[i]!=0):
      temp = state[i]
      state[i] = 0
      val = bin(i)[2:].zfill(n)
      bin_val = np.asarray(list(map(int,val)))
      bin_val_sec = bin_val.copy()
      temp2 = bin_val[m] 
      bin_val[m] = 0
      bin_val_sec[m] = 1
      new_val = int("".join(map(str,bin_val)),2)
      new_val2 = int("".join(map(str,bin_val_sec)),2)
      new_state[new_val] += temp*(1/np.sqrt(2))
      new_state[new_val2] += temp*(temp2*(-2)+1)*(1/np.sqrt(2))
  return new_state

#Implementing T Gate    
def Tgate(state,m):
  n = int(np.log2(len(state)))
  new_state = np.zeros(len(state),dtype=complex)
  for i in range(len(state)):
    if(state[i]!=0):
      temp = state[i]
      state[i] = 0
      val = bin(i)[2:].zfill(n)
      bin_val = np.asarray(list(map(int,val)))
      temp2 = bin_val[m] 
 
      new_state[i] += temp*np.exp((temp2*np.pi/4)*1j)


  return new_state

#Implementing Tdg Gate   
def TDaggate(state,m):
  n = int(np.log2(len(state)))
  new_state = np.zeros(len(state),dtype=complex)
  for i in range(len(state)):
    if(state[i]!=0):
      temp = state[i]
      state[i] = 0
      val = bin(i)[2:].zfill(n)
      bin_val = np.asarray(list(map(int,val)))
      temp2 = bin_val[m] 
      new_state[i] += temp*np.exp((-temp2*np.pi/4)*1j)

  return new_state

#Implementing Control Not Gate   
def CXgate(state,m,t):
  n = int(np.log2(len(state)))
  new_state = np.zeros(len(state),dtype=complex)
  for i in range(len(state)):
    if(state[i]!=0):
      temp = state[i]
      state[i] = 0
      val = bin(i)[2:].zfill(n)
      bin_val = np.asarray(list(map(int,val)))
      temp2 = bin_val[m]
      temp3 = bin_val[t]
      if(temp3==1):
        bin_val[m] = 1-bin_val[m]
        new_val = int("".join(map(str,bin_val)),2)
      else:
        new_val = int("".join(map(str,bin_val)),2)
      new_state[new_val] += temp
  return new_state

#Simulate function which takes a qasm string an doutputs the state vector
def simulate(qasm_string):
  global curr
  curr = -1e10
  qasm_string = qasm_string.replace('\n', '') #removing the new line literals
  qasm_string = qasm_string.split(";") #splitting on ;

  for i in range(4,len(qasm_string)):
    qasm_string[i] = qasm_string[i].replace(","," ") #replacing the ,
    temp = qasm_string[i].split(" ") # splitting on spaces

    #fetching the number of qubits in program
    if(len(temp)>0):
      for j in range(len(temp)):
        if(len(temp[j])>0 and temp[j][0]=="q"):
          l = temp[j].index('[')
          h = temp[j].index(']')
       
          curr = max(int(temp[j][l+1:h]),curr)
    
  state = np.zeros(2**(curr+1),dtype=complex)
  state[0] = 1


  # iterating over all lines and calling the gate functions
  for i in range(4,len(qasm_string)-1):
    temp = qasm_string[i].split(" ")

    if(len(temp)>0):
      op = temp[0]
      l1 = temp[1].index('[')
      h1 = temp[1].index(']')
      m1 = int(temp[1][l1+1:h1])
      if(len(temp)>2):
        l2 = temp[2].index("[")
        h2 = temp[2].index("]")
        m2 = int(temp[2][l2+1:h2])

      if(op=="x"):
        state = Xgate(state,m1)
      if(op=="h"):
        state = Hgate(state,m1)
      if(op=="t"):
        state = Tgate(state,m1)
      if(op=="tdg"):
        state = TDaggate(state,m1)
      if(op=="cx"):
        state = CXgate(state,m2,m1)
     
# rounding off since cirq output is also rounded off
  state = np.round_(state, decimals = 3)
  
  return state
 #analysis function which does various time analysis for the simulator. All pics saved in results folder
def analysis():
  global curr
  dic = {}
  dic2 = {}
  dic3 = {}
  if not os.path.exists('results'):
    os.makedirs('results')
  for i in range(3,17):
    dic[i] = 0
    dic2[i] = 0

  qasm_dir = Path(sys.argv[1])
  print("Doing Number of qubits vs time analysis")
  for qasm_file in qasm_dir.glob("**/*.qasm"):
    print("Running ",qasm_file)
    with open(qasm_file, "r") as f:
        qasm_string = f.read()

    start = time.time()
    state_vector = simulate(qasm_string)
    end = time.time()

    start2 = time.time()
    circuit = circuit_from_qasm(qasm_string)
    result = cirq.Simulator().simulate(circuit)
    statevector = list(np.around(result.state_vector(), 3))
    end2  = time.time()

  

    dic[curr+1] = end-start
    dic2[curr+1] = end2-start2
  
  


  plt.plot(list(dic.keys()),list(dic.values()))
  plt.title("Time taken vs number of qubits for my code")
  plt.xlabel("Number of qubits")
  plt.ylabel("Time Taken")
  plt.savefig("./results/time.png")  
  plt.close()

  plt.plot(list(dic2.keys()),list(dic2.values()))
  plt.title("Time taken vs number of qubits for cirq simultor")
  plt.xlabel("Number of qubits")
  plt.ylabel("Time Taken")
  plt.savefig("./results/time_cirq.png")  
  plt.close()

  print("Doing Time taken for each gate analysis")
  a1 = time.time()
  CXgate([0,0,0,1],0,1)
  a2 = time.time()
  (Hgate([1,0],0))
  a3 = time.time()
  (Tgate([0,1],0))
  a4 = time.time()
  TDaggate([0,1],0)
  a5 = time.time()
  (Xgate([1,0],0))
  a6 = time.time()
  
  dic3["X"]  = a6-a5
  dic3["H"]  = a3-a2
  dic3["T"]  = a4-a3
  dic3["Tdg"]  = a5-a4
  dic3["CX"]  = a2-a1
 
  plt.plot(dic3.keys(),dic3.values(),color = 'black')
   
  plt.title("Time taken for each Gate")
  plt.xlabel("Gate")
  plt.ylabel("Time Taken")
  plt.savefig("./results/time_gates.png") 
  plt.close()
  lis_time = []
  qasm_file = sys.argv[1]+"/miller_11.qasm"
  print("Doing Number of Gates vs Time analysis")
  with open(qasm_file, "r") as f:
      qasm_string = f.read()
      qasm_string = qasm_string.split(";")
  
      
    
      for nn in range(5,len(qasm_string)):
       
        curr = qasm_string[:nn]
        
        curr = ";".join(curr)
        b1 = time.time()
        simulate(curr)
        b2 = time.time()
        lis_time.append(b2-b1)
  
  xx = np.arange(50)
  plt.plot(xx,lis_time) 
  
  plt.xlabel("Number of Gates")
  plt.ylabel("Time Taken")
  plt.title("Time Taken vs Number of gates for miller_11.qasm")  
  plt.savefig("./results/gate_time.png")
  plt.close()    
          
def random_circuit_generate():
  
  lis = []
  lis.append('OPENQASM 2.0;\n')
  lis.append('include "qelib1.inc";\n')
  lis.append('qreg q[16];\n')
  lis.append('creg c[16];\n')
  

  qub = np.random.randint(low = 2, high = 16, size=1)[0]
  num_gates = np.random.randint(low = qub, high = 500, size=1)[0]

  if(len(sys.argv)>3):
    qub = int(sys.argv[3])
    num_gates = int(sys.argv[4])
  print("Generaing circuit with ",qub,"qubits and",num_gates," gates in total")
  print("\n")
  
  dic = {0:"h",1:"x",2:"t",3:"tdg",4:"cx"}

  for i in range(qub):
    gate = np.random.randint(low = 0, high = 5, size=1)[0]
    if(gate==4):
      qub_arr = np.arange(qub)
      qub_arr = np.delete(qub_arr, i)
      qub_num = np.random.choice(qub_arr, 1)[0]
      temp = dic[gate]+" q[{}".format(str(i))+"],q[{}".format(str(qub_num))+"]"+";\n"
      lis.append(temp)
  
    else:
      temp = dic[gate]+" q[{}]".format(str(i))+";\n" 

      lis.append(temp)
  for i in range(qub,num_gates):
    gate = np.random.randint(low = 0, high = 5, size=1)[0]
    if(gate==4):
      qub_arr = np.arange(qub)
      qub_num = np.random.choice(qub, 2, replace = False)
      temp = dic[gate]+" q[{}".format(str(qub_num[0]))+"],q[{}".format(str(qub_num[1]))+"];\n"
      lis.append(temp)
  
    else:
      qub_num = np.random.randint(low = 0, high = qub, size=1)[0]
      temp = dic[gate]+" q[{}]".format(str(qub_num))+";\n"
      lis.append(temp)
   
    
  lis_str = "".join(lis) 
  if(int(sys.argv[2])==1):
    print("Following circuit generated")
    print(lis_str)   
  state = simulate(lis_str)

  circuit = circuit_from_qasm(lis_str)
  result = cirq.Simulator().simulate(circuit)
  statevector = list(np.around(result.state_vector(), 3))
  print("Running cirq simulator and my simulator on this random circuit")
  print("Output: ",np.all(np.isclose(state, statevector)))


  return

if __name__ == '__main__':
  if(sys.argv[1]=="-gen"):
    random_circuit_generate()
  else:
    analysis()  


  