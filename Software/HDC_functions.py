import numpy as np
import random
import time

def u_gen_rand_hv(D,p=0.5):
    #use long index list, generate range and then permute them and set all numbers < p*D to 1 to get right density
    hv = [*range(D)]
    np.random.shuffle(hv)
    for i in range(len(hv)):
        if hv[i] < p*D:
          hv[i] = 1
        else:
          hv[i] = 0
    return hv


def distance(A,B):
    return sum(np.logical_xor(A,B))


def bundle_dense(block):
  if ((len(block)%2) == 0):
     block = block[0:len(block)-1]
  sums = np.sum(block, axis = 0)
  for x in range(len(sums)):
    if (sums[x] <= (len(block))/2):
      sums[x] = 0
    else:
      sums[x] = 1
  return sums

def bind_dense(A,B):
  return np.logical_xor(A,B)

def perm(A,N):
  return np.roll(A,N)


def create_item_mem_dense(N,D):
  item_mem = {}
  for n in range(N):
    item_mem[n] = u_gen_rand_hv(D)
  return item_mem


#sparse

def bundle_sparse_time_ideal(block, max_p_sparse, D):
   block = np.array(block)
   sum_ = np.sum(block,axis=0)
   output = np.array([0]*D)
   sort_index = np.argsort(-sum_) #negating array so sort is from highest to lowest
   for i in range(int(D*max_p_sparse)):
       if (sum_[sort_index[i]] == 0):
         break
       else:
          output[sort_index[i]] = 1
   return output

def bundle_sparse_time(block, threshold, D):
   sum_ = np.sum(block,axis=0)
   output = [0]*D
   for i in range(len(sum_)):
       if (sum_[i] > (threshold)): 
          output[i] = 1
   return output

def bundle_sparse_space(block):
   return np.any(block, axis=0)

def compressed_IM_sparse(N, S, L):#number of HVs, number of segments, length of segments
  item_mem = {}
  for n in range(N):
    hv = []
    for i in range(S):
      hv.append(random.randint(0,L-1))
    item_mem[n] = np.array(hv)
  return item_mem

def EM_sparse(N, S, L):
  item_mem = {}
  for n in range(N):
    hv = []
    for s in range(S):
       segment = [0]*L
       random_index = random.randint(0,L-1)
       segment[random_index] = 1
       hv.append(segment)
    item_mem[n] = np.array(hv)
  return item_mem

def binding_sparse_last(r,A,D): #r has the locations of 1's by segment in an array, A has the segments as parts of the array (2-D array!)
   rows, column_indices = np.ogrid[:A.shape[0], :A.shape[1]]
   
   # Use always a negative shift, so that column_indices are valid.
   # (could also use module operation)
   r[r < 0] += A.shape[1]
   column_indices = column_indices - r[:,np.newaxis]

   result = A[rows, column_indices] #will need to unroll to get normal vector for bundling
   return result.flatten()

def binding_sparse_segm_shift_fast(r,A,D): #r has the locations of 1's by segment in an array, A has the segments as parts of the array (2-D array!)
   rows, column_indices = np.ogrid[:A.shape[0], :A.shape[1]]

   # Use always a negative shift, so that column_indices are valid.
   # (could also use module operation)
   r += A.shape[1]
   column_indices = column_indices - r[:,np.newaxis]

   result = A[rows, column_indices] #will need to unroll to get normal vector for bundling
   return result.flatten()


def similarity_sparse(A,B,D):
   similarity = 0
   for i in range(D):
      if (A[i] == 1 & B[i] == 1):
        similarity += 1
   return similarity

def similarity_sparse_fast(A,B,D):
  A = np.array(A)
  B = np.array(B)
  return np.sum(np.logical_and(A,B))


def create_item_mem_sparse(N,D,p_sparse):
  item_mem = {}
  for n in range(N):
    item_mem[n] = u_gen_rand_hv(D,p_sparse)
  return item_mem


def binding_CDT(block_to_bind,D,random_perm_list):
  Z_hv = np.any(block_to_bind, axis=0)
  Z_hv_ORed = [0]*D
  for permut in random_perm_list:
     Z_hv_ORed = np.logical_or(Z_hv_ORed, perm(Z_hv, permut))
  return np.logical_and(Z_hv, Z_hv_ORed)
   


def random_permutation_list(K):
   result = []
   i = 0
   while (i < K):
      rand_int = random.randint(1,128)
      if rand_int not in result:
         result.append(rand_int)
         i += 1
   return result


def bundle_sparse_time_ideal_random_thinning(block, max_p_sparse_bundle, D):
  hv = np.any(block, axis=0)
  nb_1s = similarity_sparse(hv, np.array(D*[1]), D)
  while (nb_1s > max_p_sparse_bundle*D):
    i = random.randint(0, D-1)
    if (hv[i] == 1):
      hv[i] = 0
      nb_1s -= 1
  return hv


def bundle_sparse_time_ideal_CDT_thinning(block, D, random_perm_list):
  Z_hv = np.any(block, axis=0)
  Z_hv_ORed = [0]*D
  for permut in random_perm_list:
     Z_hv_ORed = np.logical_or(Z_hv_ORed, perm(Z_hv, permut))
  result = np.logical_and(Z_hv, Z_hv_ORed)
  return result

def bundle_sparse_time_ideal_laiho(block, D, nb_segments, length_segment):
  block = np.array(block)
  sum_ = np.sum(block,axis=0)
  hv_result = np.array(D*[0])
  for s in range(nb_segments):
    sort_index = np.argsort(-sum_[length_segment*s:length_segment*(s+1)])
    hv_result[sort_index[0]+s*length_segment] = 1
  return hv_result


def bundle_sparse_space_CDT(block, D, random_permut_list): #64 to bundle
  new_block = []
  for i in range(0,8):
    new_block.append(bundle_sparse_time_ideal_CDT_thinning(block[i*8:(i+1)*8], D, random_permut_list[5:6]))
  Hv = bundle_sparse_time_ideal_CDT_thinning(new_block, D, random_permut_list[0:1])
  return Hv



def bundle_sparse_time_CDT(block, D, random_permut_list):
  new_block = []
  for i in range(0,4):
    new_block.append(bundle_sparse_space_CDT(block[i*64:(i+1)*64], D, random_permut_list))

  Hv = np.any(new_block,axis=0)

  return Hv


def bundle_sparse_time_CDT2(block, D, random_permut_list):
  new_block = []
  for i in range(0,4):
    new_block.append(bundle_sparse_space_CDT(block[i*64:(i+1)*64], D, random_permut_list))
  Hv = bundle_sparse_time_ideal_CDT_thinning(new_block, D, random_permut_list[10:12])
  return Hv