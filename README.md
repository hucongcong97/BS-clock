# BS-clock
advancing epigenetic age prediction with high-resolution bisulfite sequencing data

## 1. prepare data

## 1.1 blood samples

```python
# 读取meta数据
os.chdir('G:\\Methylation\\model\\')
meta = pd.read_csv('meta.csv')
# 只保留正常样本
# meta = meta.loc[meta['disease state' == 'Normal',:]]
files = meta['SRA']
tissue = 'Blood'
for i in range(len(files)): 
    print(f"the {i}-th sample is {files[i]}")
    # 读取数据文件
    data = pd.read_csv(f'../result_new/{files[i]}.txt',sep = '\t',header=None)
    data.columns = ['Chr','Start','End','methylation-level','methylated Cs','unmethylated Cs','Strand']
    # 根据nreads筛选cpg
    data = data.loc[(data['methylated Cs']+data['unmethylated Cs'])>50,:]
    # 计算甲基化水平
    data['methylation-level'] = data['methylation-level']/100
    # 提取前四列
    data = data.iloc[:,[0,1,2,3]]
    # 根据前三列对第四列求均值
    data = data.groupby(['Chr','Start','End'])['methylation-level'].mean().reset_index()
    print(f'\tThe shape of {files[i]} is: {data.shape}')
    
    if i == 0:
        merge = data
    else:
        merge = pd.merge(merge, data, how='outer',on=['Chr','Start','End'],suffixes=(f'_x{i}',f'_y{i}'))
        # 删除重复行
        # merge = merge.drop_duplicates()
        # 删除含有NaN的行
        # merge.dropna(axis = 0,how = 'any',inplace=True)
        print(f'\t\tThe shape of merged data is: {merge.shape}')

# 设置列名
merge.columns = ['Chr','Start','End']+files.tolist()
# 删除NA值比例超过50%的行
na_count = merge.isna().sum(axis = 1)
na_ratio = na_count/(len(merge.columns)-3)
merge = merge[na_ratio <= 0.6]
print(f'**The shape of filtered merge is: {merge.shape}**')
# 保存数据
merge.to_csv(f'train_data4/merge_{tissue}_withchr.csv')
merge = merge.iloc[:,3:]

# 将merge data匹配age信息
merge_data = pd.merge(merge.T,meta[["SRA","age"]],left_index=True,right_on=["SRA"])
merge_data = merge_data.set_index("SRA")
merge_data.to_csv(f'train_data4/merge_{tissue}_withage.csv')

# 检查是否有np.nan值
merge_check = merge.T.isnull()

# 计算单个cpgsite与age的相关性
merge['corr'] = [np.corrcoef(merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),i], merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),-1])[0,1] for i in range(merge_data.shape[1]-1)]

corr_list = [0.2,0.3,0.4]
for corr in corr_list:
    data = merge[abs(merge['corr'])>corr]
    print(f'correlation threshold == {corr},the shape of data is:{data.shape}')
    # 去掉最后一列corr
    data = data.iloc[:,:-1]

    # 将NaN赋值为每个位点甲基化水平的均值
    print('Now is filling data')
    means = np.mean(data,axis=0)
    data.fillna(means, inplace=True)
    print('filling data is complete')

    data.to_csv(f'train_data4/merge_{tissue}_{corr}.csv')
```

## 1.2 Brain、Lung and Skin samples

```
# 读取meta数据
os.chdir('G:\\Methylation\\model\\')
meta = pd.read_csv('meta_tissue.csv')

# 只保留正常样本
# meta = meta.loc[meta['disease state' == 'Normal',:]]

tissue_list = ['Brain','Lung','Skin']
for tissue in tissue_list:
    print(f'Now is running {tissue}')
    meta_tissue = meta.loc[meta['tissue'] == tissue,:]
    files = meta_tissue['SRA'].reset_index(drop=True)
    
    for i in range(len(files)):
    print(f"the {i}-th sample is {files[i]}")
    # 读取数据文件
    data = pd.read_csv(f'../result_tissue/{files[i]}.txt',sep = '\t',header=None)
    data.columns = ['Chr','Start','End','methylation-level','methylated Cs','unmethylated Cs','Strand']
    # 根据nreads筛选cpg
    data = data.loc[(data['methylated Cs']+data['unmethylated Cs'])>50,:]
    # 计算甲基化水平
    data['methylation-level'] = data['methylation-level']/100
    # 提取前四列
    data = data.iloc[:,[0,1,2,3]]
    # 根据前三列对第四列求均值
    data = data.groupby(['Chr','Start','End'])['methylation-level'].mean().reset_index()
    print(f'\tThe shape of {files[i]} is: {data.shape}')

    if i == 0:
        merge = data
    else:
        merge = pd.merge(merge, data, how='outer',on=['Chr','Start','End'],suffixes=(f'_x{i}',f'_y{i}'))
        # 删除重复行
        # merge = merge.drop_duplicates()
        # 删除含有NaN的行
        # merge.dropna(axis = 0,how = 'any',inplace=True)
        print(f'\t\tThe shape of merged data is: {merge.shape}')
        
# 设置列名
merge.columns = ['Chr','Start','End']+files.tolist()
# 删除NA值比例超过50%的行
na_count = merge.isna().sum(axis = 1)
na_ratio = na_count/(len(merge.columns)-3)
merge = merge[na_ratio <= 0.6]
print(f'**The shape of filtered merge is: {merge.shape}**')
# 保存数据
merge.to_csv(f'train_data4/merge_{tissue}_withchr.csv')
merge = merge.iloc[:,3:]

# 将merge data匹配age信息
merge_data = pd.merge(merge.T,meta[["SRA","age"]],left_index=True,right_on=["SRA"])
merge_data = merge_data.set_index("SRA")
merge_data.to_csv(f'train_data4/merge_{tissue}_withage.csv')

# 检查是否有np.nan值
merge_check = merge.T.isnull()

# 计算单个cpgsite与age的相关性
merge['corr'] = [np.corrcoef(merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),i], merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),-1])[0,1] for i in range(merge_data.shape[1]-1)]

corr_list = [0.2,0.3,0.4]
for corr in corr_list:
    data = merge[abs(merge['corr'])>corr]
    print(f'correlation threshold == {corr},the shape of data is:{data.shape}')
    # 去掉最后一列corr
    data = data.iloc[:,:-1]

    # 将NaN赋值为每个位点甲基化水平的均值
    print('Now is filling data')
    means = np.mean(data,axis=0)
    data.fillna(means, inplace=True)
    print('filling data is complete')

    data.to_csv(f'train_data4/merge_{tissue}_{corr}.csv')
```

## 2. select feature

```python
def bootstrap(bootstrap_df,n):
    bootstrap_coef = []
    for i in range(n):
        #进行随机有放回抽样
        bootstrap_samples = bootstrap_df.sample(len(bootstrap_df),replace = True)
        model = ElasticNet(alpha=best_alpha, l1_ratio=best_l1_ratio, random_state=2024)
        model.fit(bootstrap_samples.iloc[:,:-1],bootstrap_samples.iloc[:,-1])
        bootstrap_coef.append(model.coef_)
    return bootstrap_coef

corr = 0.2
tissue_list = ['Blood','Brain','Lung','Skin']
for tissue in tissue_list:
    if tissue == 'Blood':
        meta = pd.read_csv('meta.csv')
    else:
        df = pd.read_csv('meta_tissue.csv')
        meta = df.loc[df['tissue'] == tissue,:]
    data = pd.read_csv(f'train_data4/merge_{tissue}_{corr}.csv',index_col=0)
    data = data.T
    # 填充缺失值
    means = np.mean(data,axis=0)
    data.fillna(means, inplace=True)

    # 利用10-fold交叉验证求出最优参数
    X = np.array(data)
    bootstrap_data = pd.merge(data,meta[["SRA","age"]],left_index=True,right_on=["SRA"])
    bootstrap_data = bootstrap_data.set_index("SRA")
    y = bootstrap_data['age']
    
    # 定义参数搜索范围
    alpha_range = [0.1,1,10,50]
    l1_ratio_range = np.arange(0, 1, 0.1)

    best_alpha = None
    best_l1_ratio = None
    best_mean_mse = float('inf')
    # Best alpha: 0.1, Best l1_ratio: 0.4, Best mean MSE: 97.13731881603614

    # 交叉验证选参数
    print('Now is running 10-fold cross validation')
    t1 = time.time()
    for alpha in alpha_range:
        for l1_ratio in l1_ratio_range:
            # 定义 K 折交叉验证
            k_folds = KFold(n_splits=10, shuffle=True, random_state=2024)

            # 存储每次交叉验证的均方误差
            mean_mse_list = []

            # 使用 K 折交叉验证训练模型并进行评估
            for train_index, test_index in k_folds.split(X):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]

                # 创建并训练模型
                model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=2024)
                model.fit(X_train, y_train)

                # 在测试集上进行预测
                y_pred = model.predict(X_test)

                # 计算均方误差
                mse = mean_squared_error(y_test, y_pred)
                mean_mse_list.append(mse)

            # 计算交叉验证均方误差的平均值
            mean_mse = np.mean(mean_mse_list)

            # 更新最佳参数和最小均方误差
            if mean_mse < best_mean_mse:
                best_mean_mse = mean_mse
                best_alpha = alpha
                best_l1_ratio = l1_ratio

    #           print(f"alpha = {round(alpha, 2)}, l1_ratio = {round(l1_ratio, 3)}, mean MSE = {round(mean_mse, 3)}")

    # 输出最佳参数
    print(f"Best alpha: {best_alpha}, Best l1_ratio: {best_l1_ratio}, Best mean MSE: {best_mean_mse}")
    t2 = time.time()
    print(f'10-fold cross validation is complete,takes {str(t2-t1)}')

    # bootsrtap得到500个模型，并得到弹性网回归的系数
    print('Now is bootstrap')
    t3 = time.time()
    coef = bootstrap(bootstrap_data,500)
    coef = pd.DataFrame(coef,columns=list(data))
    t4 = time.time()
    print(f'bootstrap is complete,takes {str(t4-t3)}')

    # 删除在少于一半的模型（250）中选择的特征
    zero_percent = (coef == 0).sum(axis = 0) / 500
    columns_to_drop = zero_percent[zero_percent > 0.5].index

    # 输出筛选后500个模型中重要特征的系数
    coef_drop = coef.drop(columns_to_drop,axis = 1)

    # 确定最后模型所要使用的数据
    model_data = bootstrap_data.drop(columns_to_drop ,axis=1)
    print(f'The shape of model_data of {tissue} is: {model_data.shape}')
    model_data.to_csv(f'train_data4/bootstrap_{tissue}_{corr}.csv') 
```

## 3. train model

### 3.1 without disease features

```python
# 定义参数搜索范围
alpha_range = [0.1,1,10,50]
l1_ratio_range = np.arange(0, 1, 0.05)

corr = 0.2
tissue_list = ['Blood','Brain','Lung','Skin']
for tissue in tissue_list:
    model_data = pd.read_csv(f'train_data4/bootstrap_{tissue}_{corr}.csv',index_col=0)

    loo = LeaveOneOut()
    # 存储预测的年龄
    age_pred_list = []
    # 使用留一法分割数据集
    print('Now is running LOOCV')
    t1 = time.time()
    for train,test in loo.split(model_data):
        train_x = model_data.iloc[train].iloc[:,:-1]
        train_y = model_data.iloc[train].iloc[:,-1]
        test_x = model_data.iloc[test].iloc[:,:-1]
        test_y = model_data.iloc[test].iloc[:,-1]

        best_alpha = None
        best_l1_ratio = None
        best_mean_mse = float('inf')

        # K 折交叉验证选择最优参数
        for alpha in alpha_range:
            for l1_ratio in l1_ratio_range:
                # 定义 K 折交叉验证
                k_folds = KFold(n_splits=5, shuffle=True, random_state=2024)

                # 存储每次交叉验证的均方误差
                mean_mse_list = []
                # 使用 K 折交叉验证训练模型并进行评估
                for train_index, test_index in k_folds.split(train_x):
                    X_train, X_test = train_x.iloc[train_index], train_x.iloc[test_index]
                    y_train, y_test = train_y.iloc[train_index], train_y.iloc[test_index]

                    # 创建并训练模型
                    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=2024)
                    model.fit(X_train, y_train)

                    # 在测试集上进行预测
                    y_pred = model.predict(X_test)

                    # 计算均方误差
                    mse = mean_squared_error(y_test, y_pred)
                    mean_mse_list.append(mse)

                # 计算交叉验证均方误差的平均值
                mean_mse = np.mean(mean_mse_list)

                # 更新最佳参数和最小均方误差
                if mean_mse < best_mean_mse:
                    best_mean_mse = mean_mse
                    best_alpha = alpha
                    best_l1_ratio = l1_ratio

#                     print(f"alpha = {round(alpha, 2)}, l1_ratio = {round(l1_ratio, 3)}, mean MSE = {round(mean_mse, 3)}")

        # 输出最佳参数
        print(f"Best alpha: {best_alpha}, Best l1_ratio: {best_l1_ratio}, Best mean MSE: {best_mean_mse}")

        # 留一法预测年龄
    #     best_alpha = 0.1
    #     best_l1_ratio = 0.4

        age_model = ElasticNet(alpha=best_alpha, l1_ratio=best_l1_ratio, random_state=2024)
        age_model.fit(train_x,train_y)
        age_pred = age_model.predict(test_x)
        age_pred_list.append(age_pred[0])

    t2 = time.time()
    print(f'LOOCV step takes {str(t2-t1)}')
    df = pd.DataFrame(age_pred_list)
    df.to_csv(f'train_data4/age_{tissue}_{corr}.csv')
    model_corr = np.corrcoef(age_pred_list,model_data["age"])[0, 1]
    print(f'model_corr of {tissue} is:{model_corr}')
```

### 3.2 with disease features

```python
# add disease feature
corr = 0.2
tissue_list = ['Blood','Brain','Lung','Skin']
for tissue in tissue_list:
    print(tissue)
    if tissue == 'Blood':
        meta = pd.read_csv('train_data4/meta.csv')
    else:
        df = pd.read_csv('train_data4/meta_tissue.csv')
        meta = df.loc[df['tissue'] == tissue,:]

    model_data = pd.read_csv(f'train_data4/bootstrap_{tissue}_{corr}.csv',index_col=0)
    # 根据meta将样本重排
    model_data = model_data.loc[meta['SRA']]
    # 根据疾病信息增加列
    for i in set(meta['disease state']):
        model_data[i] = (meta['disease state']==i).astype(int).tolist()
    # 将age列调整为最后一列
    model_data = model_data[[col for col in model_data.columns if col != 'age'] + ['age']]
    # 保存数据
    model_data.to_csv(f'train_data4/bootstrap_{tissue}_{corr}_add.csv') 
```

```
# train model
# 定义参数搜索范围
alpha_range = [0.1,1,10,50]
l1_ratio_range = np.arange(0, 1, 0.05)

corr = 0.2
tissue_list = ['Blood','Brain','Lung','Skin']
for tissue in tissue_list:
    model_data = pd.read_csv(f'train_data4/bootstrap_{tissue}_{corr}_add.csv',index_col=0)

    loo = LeaveOneOut()
    # 存储预测的年龄
    age_pred_list = []
    # 使用留一法分割数据集
    print('Now is running LOOCV')
    t1 = time.time()
    for train,test in loo.split(model_data):
        train_x = model_data.iloc[train].iloc[:,:-1]
        train_y = model_data.iloc[train].iloc[:,-1]
        test_x = model_data.iloc[test].iloc[:,:-1]
        test_y = model_data.iloc[test].iloc[:,-1]

        best_alpha = None
        best_l1_ratio = None
        best_mean_mse = float('inf')

        # K 折交叉验证选择最优参数
        for alpha in alpha_range:
            for l1_ratio in l1_ratio_range:
                # 定义 K 折交叉验证
                k_folds = KFold(n_splits=5, shuffle=True, random_state=2024)

                # 存储每次交叉验证的均方误差
                mean_mse_list = []
                # 使用 K 折交叉验证训练模型并进行评估
                for train_index, test_index in k_folds.split(train_x):
                    X_train, X_test = train_x.iloc[train_index], train_x.iloc[test_index]
                    y_train, y_test = train_y.iloc[train_index], train_y.iloc[test_index]

                    # 创建并训练模型
                    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=2024)
                    model.fit(X_train, y_train)

                    # 在测试集上进行预测
                    y_pred = model.predict(X_test)

                    # 计算均方误差
                    mse = mean_squared_error(y_test, y_pred)
                    mean_mse_list.append(mse)

                # 计算交叉验证均方误差的平均值
                mean_mse = np.mean(mean_mse_list)

                # 更新最佳参数和最小均方误差
                if mean_mse < best_mean_mse:
                    best_mean_mse = mean_mse
                    best_alpha = alpha
                    best_l1_ratio = l1_ratio

#                     print(f"alpha = {round(alpha, 2)}, l1_ratio = {round(l1_ratio, 3)}, mean MSE = {round(mean_mse, 3)}")

        # 输出最佳参数
        print(f"Best alpha: {best_alpha}, Best l1_ratio: {best_l1_ratio}, Best mean MSE: {best_mean_mse}")

        # 留一法预测年龄
    #     best_alpha = 0.1
    #     best_l1_ratio = 0.4

        age_model = ElasticNet(alpha=best_alpha, l1_ratio=best_l1_ratio, random_state=2024)
        age_model.fit(train_x,train_y)
        age_pred = age_model.predict(test_x)
        age_pred_list.append(age_pred[0])

    t2 = time.time()
    print(f'LOOCV step takes {str(t2-t1)}')
    df = pd.DataFrame(age_pred_list)
    df.to_csv(f'train_data4/age_{tissue}_{corr}_add.csv')
    model_corr = np.corrcoef(age_pred_list,model_data["age"])[0, 1]
    print(f'model_corr of {tissue} is:{model_corr}')
```

