import numpy as np
from numpy.linalg import solve
import pandas as pd
import time
import sys
import csv
import math
import datetime
import itertools
import psutil
import random
import os
import networkx as nx


def getTransitionProbability(N, R = 1):
    pr = np.zeros((R * N, R * N))
    for r in range(R):
        class_number = 0
        p = np.random.rand(N, N)
        for i, val in enumerate(np.sum(p, axis=1)):
            p[i] /= val
        for i in range(N):
            for j in range(N):
                pr[r * N+i,r * N+j] = p[i,j]
    return pr

def getEquivalence4(th, roop, p):
  reachable = np.zeros((len(p), len(p))) #全て0のpと同じ大きさの行列を用意
  equivalence = [[] for i in range(len(p))] #同値類を最大数用意
  list_number = 0

  for ix in range(roop): #n乗まで実施
    pn = np.linalg.matrix_power(p.copy(), ix+1) #累乗 
    for i in range(len(pn)):
      for j in range(len(pn)):
        if(pn[i][j] > th): #推移確率が閾値より大きいかチェック
          reachable[i,j] = 1

#  print('到達行列 = \n{0}'.format(reachable))
#  communicate = np.zeros((len(p), len(p))) #全て0のpと同じ大きさの行列を用意
#  for i in range(len(reachable)):
#    for j in range(len(reachable)):
#      if(reachable[i][j] == 1 and reachable[j][i] == 1):
    
  for i in range(len(pn)):
    for j in range(i+1, len(pn)):
      if(reachable[i][j] == 1 and reachable[j][i] == 1):
#        print('reachable {0} <-> {1}'.format(i, j))
        find = 0 #iでfind -> 1, jでfind -> 2
        idx = len(pn)
        for k in range(len(pn)):
          if i in equivalence[k]: #iが同値類k番目に存在している
            find = 1 #iは同値類k番目に存在
            idx = k
#            print('{0} find in {1}'.format(i, equivalence[k]))
            break
          elif j in equivalence[k]: 
            find = 2
            idx = k
#            print('{0} find in {1}'.format(j, equivalence[k]))
            break
        if find == 1:
          if j not in equivalence[idx]: #jは同値類kには存在しない (他のリストにもないことを確認する!!!他のリストにあった場合は移動がいい？ -> communicateがずれないように最後に集合で演算する)
            equivalence[idx].append(j) #jを追加
#            print('{0}に{1}を追加'.format(equivalence[idx], j))
#            print('リスト全体 {0}'.format(equivalence))     
            #break   
        elif find == 2:
          if i not in equivalence[idx]: #(他のリストにもないことを確認する!!!)
            equivalence[idx].append(i)
#            print('{0}に{1}を追加'.format(equivalence[idx], i))
#            print('リスト全体 {0}'.format(equivalence))         
            #break
        elif(find == 0): #新規に追加
          equivalence[list_number].append(i)
#          print('{0}に{1}を追加'.format(equivalence[list_number], i))
#          print('リスト全体 {0}'.format(equivalence))
          if(i != j):
            equivalence[list_number].append(j)
#            print('{0}に{1}を追加'.format(equivalence[list_number], j))
#            print('リスト全体 {0}'.format(equivalence))
          list_number += 1

  #4. Communicateにならないノードを登録
  for i in range(len(pn)):
    find = 0
    for j in range(len(pn)):
      if i in equivalence[j]:
        find = 1
        break
    if find == 0:
      equivalence[list_number].append(i)
#      print('Non Communication node {0} add'.format(i))
#      print('リスト全体 {0}'.format(equivalence))
      list_number += 1

  #5. エルゴード性の確認(class数が1のとき)
  class_number = 0
  for i in range(len(pn)):
    if len(equivalence[i]) > 0:
#      print('クラスの長さ : {0}, {1}'.format(len(equivalence[i]), equivalence[i]))
      class_number += 1

  for i in range(class_number):
    for j in range(i+1, class_number):
      s1 = set(equivalence[i])
      s2 = set(equivalence[j])
      if s1 & s2 :
#        print('重複 {0} & {1}'.format(i, j))
#        print('重複ノード {}'.format(s1 & s2))
        equivalence[i] = equivalence[i] + equivalence[j]
        equivalence[j].clear()

  #print('修正クラス数 {0}'.format(modify_class_number))
  for i in range(class_number):
  #  print(equivalence[i])
    equivalence[i] =list(set(equivalence[i]))

  #再度クラス数チェック
  class_number = 0
  for i in range(len(pn)):
    if len(equivalence[i]) > 0:
      class_number += 1
  
  #再度リストを構成
  modify_equivalence = [[] for i in range(class_number)] #同値類をクラス数用意
  l_index = 0
  for i in range(len(pn)):
    if len(equivalence[i]) > 0:
      modify_equivalence[l_index] = equivalence[i]
      l_index += 1

  return modify_equivalence, class_number, reachable

def getDegreeMax(equivalence, reachable):
  degree_max_value = [] #クラスの中で最大次数のものをクラス0から順番に入れる
  degree_max_index = []
  reachable_sum = np.sum(reachable, axis=0) + np.sum(reachable, axis=1) #到達行列から行和と列和の和
  for i in range(len(equivalence)): #クラス数だけ回す
    deg_max = -1 #最大値の初期値
    deg_max_index = -1 #最大値を持つノード番号
    for j in equivalence[i]: #ノード番号を一つずつ取り出す
      if deg_max < reachable_sum[j]:
        deg_max = reachable_sum[j]
        deg_index = j
    degree_max_value.append(deg_max)
    degree_max_index.append(deg_index)
  return degree_max_index, degree_max_value

def getNormalize(p, degree_max_index, p_plus = 0.01):
  for i in degree_max_index:
    for j in degree_max_index:
      if i != j: #(i,j),(j,i)両方に足す
        p[i][j] += p_plus 
  #推移確率に直す
  for i, val in enumerate(np.sum(p, axis=1)):
    p[i] /= val
  return p

def getAdjacency(reachable, equi, delta):
  adjacency = [[0 for i in range(len(equi))] for j in range(len(equi))]
  for i in range(len(equi)):
    for j in range(len(equi)):
      if reachable[equi[i]][equi[j]] >= delta:
        adjacency[i][j] = 1
  #print(adjacency)
  return adjacency

def getProximityCentrality(equivalence, reachable):
  degree_max_value = [] #クラスの中で最大次数のものをクラス0から順番に入れる
  degree_max_index = []
  for i in range(len(equivalence)): #クラス数だけ回す
    if len(equivalence[i]) > 1:
      edges = []
      for e1 in equivalence[i]:
        for e2 in equivalence[i]:
          if reachable[e1][e2] == 1:
            edges.append((e1,e2))
      G = nx.DiGraph(edges)
      centrality = nx.closeness_centrality(G)
      deg_max = max(centrality.values())
      deg_index = max(centrality, key=centrality.get)
      #print(max(centrality, key=centrality.get),max(centrality.values()), centrality,deg_index,deg_max)
      #print(equivalence[i])
    else:
      deg_max = 0
      deg_index = equivalence[i][0]
    degree_max_value.append(deg_max)
    degree_max_index.append(deg_index)
  return degree_max_index, degree_max_value
  #1隣接行列を作成, 変数定義(次数変数)

def getMediaCentricity(equivalence, reachable):#媒介中心性, delta
  degree_max_value = [] #クラスの中で最大次数のものをクラス0から順番に入れる
  degree_max_index = []
  
  #reachable_sum = np.sum(reachable, axis=0) + np.sum(reachable, axis=1) #到達行列から行和と列和の和
  for i in range(len(equivalence)): #クラス数だけ回す
    if len(equivalence[i]) > 1:
      edges = []
      for e1 in equivalence[i]:
        for e2 in equivalence[i]:
          if reachable[e1][e2] == 1:
            edges.append((e1,e2))
      G = nx.DiGraph(edges)
      centrality = nx.betweenness_centrality(G)
      deg_max = max(centrality.values())
      deg_index = max(centrality, key=centrality.get)
      #print(max(centrality, key=centrality.get),max(centrality.values()), centrality,deg_index,deg_max)
      #print(equivalence[i])
    else:
      deg_max = 0
      deg_index = equivalence[i][0]
    degree_max_value.append(deg_max)
    degree_max_index.append(deg_index)
  return degree_max_index, degree_max_value

def getEigenvectorCentricity(equivalence, reachable):#固有ベクトル中心性, delta
  degree_max_value = [] #クラスの中で最大次数のものをクラス0から順番に入れる
  degree_max_index = []
  
  #reachable_sum = np.sum(reachable, axis=0) + np.sum(reachable, axis=1) #到達行列から行和と列和の和
  for i in range(len(equivalence)): #クラス数だけ回す
    if len(equivalence[i]) > 1:
      edges = []
      for e1 in equivalence[i]:
        for e2 in equivalence[i]:
          if reachable[e1][e2] == 1:
            edges.append((e1,e2))
      G = nx.DiGraph(edges)
      centrality = nx.eigenvector_centrality(G)
      deg_max = max(centrality.values())
      deg_index = max(centrality, key=centrality.get)
      #print(max(centrality, key=centrality.get),max(centrality.values()), centrality,deg_index,deg_max)
      #print(equivalence[i])
    else:
      deg_max = 0
      deg_index = equivalence[i][0]
    degree_max_value.append(deg_max)
    degree_max_index.append(deg_index)
  return degree_max_index, degree_max_value


def getPageRank(equivalence, reachable):#PageRank, delta
  degree_max_value = [] #クラスの中で最大次数のものをクラス0から順番に入れる
  degree_max_index = []
  
  #reachable_sum = np.sum(reachable, axis=0) + np.sum(reachable, axis=1) #到達行列から行和と列和の和
  for i in range(len(equivalence)): #クラス数だけ回す
    if len(equivalence[i]) > 1:
      edges = []
      for e1 in equivalence[i]:
        for e2 in equivalence[i]:
          if reachable[e1][e2] == 1:
            edges.append((e1,e2))
      G = nx.DiGraph(edges)
      centrality = nx.pagerank(G)
      deg_max = max(centrality.values())
      deg_index = max(centrality, key=centrality.get)
      #print(max(centrality, key=centrality.get),max(centrality.values()), centrality,deg_index,deg_max)
      #print(equivalence[i])
    else:
      deg_max = 0
      deg_index = equivalence[i][0]
    degree_max_value.append(deg_max)
    degree_max_index.append(deg_index)
  return degree_max_index, degree_max_value


def getInformationCentricity(equivalence, reachable):#情報中心性 無向グラフのみ, delta
  degree_max_value = [] #クラスの中で最大次数のものをクラス0から順番に入れる
  degree_max_index = []
  
  #reachable_sum = np.sum(reachable, axis=0) + np.sum(reachable, axis=1) #到達行列から行和と列和の和
  for i in range(len(equivalence)): #クラス数だけ回す
    if len(equivalence[i]) > 1:
      edges = []
      for e1 in equivalence[i]:
        for e2 in equivalence[i]:
          if reachable[e1][e2] == 1:
            edges.append((e1,e2))
      G = nx.Graph(edges)
      centrality = nx.information_centrality(G)
      deg_max = max(centrality.values())
      deg_index = max(centrality, key=centrality.get)
      #print(max(centrality, key=centrality.get),max(centrality.values()), centrality,deg_index,deg_max)
      #print(equivalence[i])
    else:
      deg_max = 0
      deg_index = equivalence[i][0]
    degree_max_value.append(deg_max)
    degree_max_index.append(deg_index)
  return degree_max_index, degree_max_value


def savedata(N, delta, p_plus, roop_num, dpath, name, class_number_list, class_each_number_list, degree_max_value_list, pr, p):
    
    path = './csv/{}/'.format(dpath)
    file = str(d)+'_'+str(name)+'_'+str(N)+'_delta_'+str(delta)+'_p_plus_'+str(p_plus)+'_roop_num_'+str(roop_num)+'_class_number_list.csv'
    np.savetxt(path + file, np.array(class_number_list), delimiter=",", fmt="%d")

    file = str(d)+'_'+str(name)+'_'+str(N)+'_delta_'+str(delta)+'_p_plus_'+str(p_plus)+'_roop_num_'+str(roop_num)+'_class_each_number_list.csv'
    with open(path + file, 'w') as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerows(class_each_number_list)

    file = str(d)+'_'+str(name)+'_'+str(N)+'_delta_'+str(delta)+'_p_plus_'+str(p_plus)+'_roop_num_'+str(roop_num)+'_degree_max_value_list.csv'
    with open(path + file, 'w') as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerows(degree_max_value_list)

    file = str(d)+'_'+str(name)+'_'+str(N)+'_delta_'+str(delta)+'_p_plus_'+str(p_plus)+'_roop_num_'+str(roop_num)+'_initial_p.csv'
    np.savetxt(path + file, pr, delimiter=",", fmt="%f")

    file = str(d)+'_'+str(name)+'_'+str(N)+'_delta_'+str(delta)+'_p_plus_'+str(p_plus)+'_roop_num_'+str(roop_num)+'_last_p.csv'
    np.savetxt(path + file, p, delimiter=",", fmt="%f")

def append_summarydata(summarydata, sroop, degree_max_index, degree_max_value, equivalencelist, degree_max_value_list):
    summarydata[0].append(sroop)
    summarydata[1].append(len(degree_max_index))
    summarydata[2].append(degree_max_index[0] if len(degree_max_index) == 1 else degree_max_index)
    summarydata[3].append(degree_max_value[0] if len(degree_max_index) == 1 else degree_max_value)
    summarydata[4].append(equivalencelist)
    summarydata[5].append(degree_max_value_list)
    return summarydata

def calculate_roop(pr, roop_num, delta, N, d, summarydata, calflag):
  #calflagname = ['■■次数中心性■■', '■■近接中心性■■', '■■媒介中心性■■', '■■固有ベクトル中心性■■', '■■PageRank■■']
  calflagname = ['DegreeMax', 'ProximityCentrality', 'MediaCentricity', 'EigenvectorCentricity', 'PageRank']
  name = calflagname[calflag]
  print(name)
  roop = 1
  sroop = 1
  equivalencelist = []
  p = pr.copy()
  class_number_list = [] #各回のクラス数を保存
  degree_max_value_list = [] #各回のクラス毎の最大次数を保存
  class_each_number_list = [[] for i in range(roop_num)] #各回のそれぞれのクラスの要素数
  for i in range(roop_num):
    print('{0}回目'.format(i))
    equivalence, class_number, reachable = getEquivalence4(delta, roop, p)
    equivalencelist.append(equivalence)
    print('クラス数 : {0}'.format(class_number))
    class_number_list.append(class_number)
    if class_number != N:
      roopflag = 0
    else:
      roopflag += 1
      if roopflag == 50:
        sroop = i * -1
        break
    for j in range(class_number):
      #print(equivalence[j])
      class_each_number_list[i].append(len(equivalence[j]))
    if calflag == 0:
        degree_max_index, degree_max_value = getDegreeMax(equivalence, reachable)
    elif calflag == 1:
        degree_max_index, degree_max_value = getProximityCentrality(equivalence, reachable)
    elif calflag == 2:
        degree_max_index, degree_max_value = getMediaCentricity(equivalence, reachable)
    elif calflag == 3:
        try :
            degree_max_index, degree_max_value = getEigenvectorCentricity(equivalence, reachable)
        except:
            sroop = i * -1
            break
    elif calflag == 4:
        degree_max_index, degree_max_value = getPageRank(equivalence, reachable)
    else :
        print('error')
    #print('各クラスの最大次数を持つノード番号 : {0}'.format(degree_max_index))
    #print('各クラスの最大次数 : {0}'.format(degree_max_value))
    degree_max_value_list.append(degree_max_value)
    p = getNormalize(p, degree_max_index)
    sroop = i
    if class_number == 1:
      break
    #print('\n')
  #print(class_number_list)
  #print(class_each_number_list)
  #print(degree_max_value_list)
  #name = 'DegreeMax'
  savedata(N, delta, p_plus, roop_num, d, name, class_number_list, class_each_number_list, degree_max_value_list, pr, p)
  summarydata = append_summarydata(summarydata, sroop, degree_max_index, degree_max_value, equivalencelist, degree_max_value_list)
  return summarydata

N = int(sys.argv[1])#33
delta = 1 / N * float(sys.argv[2])#1.5
#delta = 0.045
p_plus = float(sys.argv[3])#0.001
roop_num = int(sys.argv[4])#20
pr = getTransitionProbability(N)

summarydata = [[],[],[],[],[],[]]

#データの保存
t_delta = datetime.timedelta(hours=9)
JST = datetime.timezone(t_delta, 'JST')
now = datetime.datetime.now(JST)
d = now.strftime('%Y%m%d%H%M%S')
d = '{}_{}_{}_{}_{}'.format(N, delta, p_plus, roop_num, d)
os.mkdir('./csv/{}'.format(d))
  
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 0)#DegreeMax
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 1)#ProximityCentrality
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 2)#MediaCentricity'
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 3)#EigenvectorCentricity
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 4)#PageRank

#print(summarydata)
with open('./csv/{}/summary.csv'.format(d), 'w') as file:
    writer = csv.writer(file, lineterminator='\n')
    writer.writerows(summarydata)
#print('./csv/{}'.format(d))


summarydatalow = [sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]] + list(itertools.chain.from_iterable(summarydata))#summarydata[0] + summarydata[1] + summarydata[2] + summarydata[3]+summarydata[4]+summarydata[5]
#print(summarydatalow)
with open('./summary.csv'.format(d), 'a') as file:
    writer = csv.writer(file, lineterminator='\n')
    writer.writerow(summarydatalow)

print(summarydata[0])
#python3 DMC_Connection.py 33 1.5 0.001 20