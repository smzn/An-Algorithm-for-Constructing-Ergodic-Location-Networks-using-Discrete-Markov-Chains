Program to connect networks using centrality so that the number of equivalence classes is 1 from the initial transition probability (random occurrence)

Execution example: python3 DMC_Connection.py 33 1.5 0.001 20
 　　　　python3 DMC_Connection.py (number of nodes) (delta threshold coefficient} (ε) (roop number)

Details of the results are stored in the csv folder

For lines 405~409
 　The last value determines the centrality to be used
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 0)#DegreeMax
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 1)#ProximityCentrality
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 2)#MediaCentricity'
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 3)#EigenvectorCentricity
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 4)#PageRank





初期推移確率（ランダム発生）から同値類数が1になるように、中心性を用いてネットワークを接続していくプログラム

実行例：python3 DMC_Connection.py 33 1.5 0.001 20
 　　　　python3 DMC_Connection.py (ノード数) (delta閾値係数) (ε) (roop数)

結果の詳細はcsvフォルダ内に保存されます

プログラム内405~409のプログラムについて
 　最後の値によって利用する中心性を決めています
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 0)#DegreeMax
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 1)#ProximityCentrality
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 2)#MediaCentricity'
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 3)#EigenvectorCentricity
summarydata = calculate_roop(pr, roop_num, delta, N, d, summarydata, 4)#PageRank


