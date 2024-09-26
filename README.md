# BS-clock
advancing epigenetic age prediction with high-resolution bisulfite sequencing data

## 1. prepare data

## 1.1 blood samples

```python
# Set the working directory
os.chdir('G:\\Methylation\\model\\')
meta = pd.read_csv('meta.csv')
# Only keep normal samples
# meta = meta.loc[meta['disease state' == 'Normal',:]]
files = meta['SRA']
tissue = 'Blood'
for i in range(len(files)): 
    print(f"the {i}-th sample is {files[i]}")
    # read data
    data = pd.read_csv(f'../result_new/{files[i]}.txt',sep = '\t',header=None)
    data.columns = ['Chr','Start','End','methylation-level','methylated Cs','unmethylated Cs','Strand']
    # Filter CpGs based on nreads
    data = data.loc[(data['methylated Cs']+data['unmethylated Cs'])>50,:]
    # Calculate mean methylation level
    data['methylation-level'] = data['methylation-level']/100
    # Extract the first four columns
    data = data.iloc[:,[0,1,2,3]]
    # Calculate the average of the fourth column based on the first three columns
    data = data.groupby(['Chr','Start','End'])['methylation-level'].mean().reset_index()
    print(f'\tThe shape of {files[i]} is: {data.shape}')
    
    if i == 0:
        merge = data
    else:
        merge = pd.merge(merge, data, how='outer',on=['Chr','Start','End'],suffixes=(f'_x{i}',f'_y{i}'))
        # Delete duplicate rows
        # merge = merge.drop_duplicates()
        # Delete rows containing NaN
        # merge.dropna(axis = 0,how = 'any',inplace=True)
        print(f'\t\tThe shape of merged data is: {merge.shape}')

# Set column names
merge.columns = ['Chr','Start','End']+files.tolist()
# Delete rows with NA value ratio exceeding 50%
na_count = merge.isna().sum(axis = 1)
na_ratio = na_count/(len(merge.columns)-3)
merge = merge[na_ratio <= 0.6]
print(f'**The shape of filtered merge is: {merge.shape}**')
# Save data
merge.to_csv(f'train_data4/merge_{tissue}_withchr.csv')
merge = merge.iloc[:,3:]

# Match merge data to age information
merge_data = pd.merge(merge.T,meta[["SRA","age"]],left_index=True,right_on=["SRA"])
merge_data = merge_data.set_index("SRA")
merge_data.to_csv(f'train_data4/merge_{tissue}_withage.csv')

# Check for np.nan values
merge_check = merge.T.isnull()

# Calculate the correlation between single CpG site and age
merge['corr'] = [np.corrcoef(merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),i], merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),-1])[0,1] for i in range(merge_data.shape[1]-1)]

corr_list = [0.2,0.3,0.4]
for corr in corr_list:
    data = merge[abs(merge['corr'])>corr]
    print(f'correlation threshold == {corr},the shape of data is:{data.shape}')
    # Remove the last column corr
    data = data.iloc[:,:-1]

    # Replace NaN with the mean of methylation level at each site
    print('Now is filling data')
    means = np.mean(data,axis=0)
    data.fillna(means, inplace=True)
    print('filling data is complete')

    data.to_csv(f'train_data4/merge_{tissue}_{corr}.csv')
```

## 1.2 Brain、Lung and Skin samples

```
# Set the working directory
os.chdir('G:\\Methylation\\model\\')
meta = pd.read_csv('meta_tissue.csv')

# Only keep normal samples
# meta = meta.loc[meta['disease state' == 'Normal',:]]

tissue_list = ['Brain','Lung','Skin']
for tissue in tissue_list:
    print(f'Now is running {tissue}')
    meta_tissue = meta.loc[meta['tissue'] == tissue,:]
    files = meta_tissue['SRA'].reset_index(drop=True)
    
    for i in range(len(files)):
    print(f"the {i}-th sample is {files[i]}")
    # read data
    data = pd.read_csv(f'../result_tissue/{files[i]}.txt',sep = '\t',header=None)
    data.columns = ['Chr','Start','End','methylation-level','methylated Cs','unmethylated Cs','Strand']
    # Filter CpGs based on nreads
    data = data.loc[(data['methylated Cs']+data['unmethylated Cs'])>50,:]
    # Calculate mean methylation level
    data['methylation-level'] = data['methylation-level']/100
    # Extract the first four columns
    data = data.iloc[:,[0,1,2,3]]
    # Calculate the average of the fourth column based on the first three columns
    data = data.groupby(['Chr','Start','End'])['methylation-level'].mean().reset_index()
    print(f'\tThe shape of {files[i]} is: {data.shape}')

    if i == 0:
        merge = data
    else:
        merge = pd.merge(merge, data, how='outer',on=['Chr','Start','End'],suffixes=(f'_x{i}',f'_y{i}'))
        # Set column names
        # merge = merge.drop_duplicates()
        # Delete rows containing NaN
        # merge.dropna(axis = 0,how = 'any',inplace=True)
        print(f'\t\tThe shape of merged data is: {merge.shape}')
        
# Set column names
merge.columns = ['Chr','Start','End']+files.tolist()
# Delete rows with NA value ratio exceeding 50%
na_count = merge.isna().sum(axis = 1)
na_ratio = na_count/(len(merge.columns)-3)
merge = merge[na_ratio <= 0.6]
print(f'**The shape of filtered merge is: {merge.shape}**')
# Save data
merge.to_csv(f'train_data4/merge_{tissue}_withchr.csv')
merge = merge.iloc[:,3:]

# Match merge data to age information
merge_data = pd.merge(merge.T,meta[["SRA","age"]],left_index=True,right_on=["SRA"])
merge_data = merge_data.set_index("SRA")
merge_data.to_csv(f'train_data4/merge_{tissue}_withage.csv')

# Check for np.nan values
merge_check = merge.T.isnull()

# Calculate the correlation between single CpG site and age
merge['corr'] = [np.corrcoef(merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),i], merge_data.iloc[(~merge_check.iloc[:,i]).tolist(),-1])[0,1] for i in range(merge_data.shape[1]-1)]

corr_list = [0.2,0.3,0.4]
for corr in corr_list:
    data = merge[abs(merge['corr'])>corr]
    print(f'correlation threshold == {corr},the shape of data is:{data.shape}')
    # Remove the last column corr
    data = data.iloc[:,:-1]

    # Replace NaN with the mean of methylation level at each site
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
        # Perform random sampling with replacement
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
    # Fill missing values
    means = np.mean(data,axis=0)
    data.fillna(means, inplace=True)

    # Find optimal parameters using 10-fold cross validation
    X = np.array(data)
    bootstrap_data = pd.merge(data,meta[["SRA","age"]],left_index=True,right_on=["SRA"])
    bootstrap_data = bootstrap_data.set_index("SRA")
    y = bootstrap_data['age']
    
    # Define parameter search range
    alpha_range = [0.1,1,10,50]
    l1_ratio_range = np.arange(0, 1, 0.1)

    best_alpha = None
    best_l1_ratio = None
    best_mean_mse = float('inf')
    # Best alpha: 0.1, Best l1_ratio: 0.4, Best mean MSE: 97.13731881603614

    # Cross validation parameter selection
    print('Now is running 10-fold cross validation')
    t1 = time.time()
    for alpha in alpha_range:
        for l1_ratio in l1_ratio_range:
            # Define K-fold cross validation
            k_folds = KFold(n_splits=10, shuffle=True, random_state=2024)

            # Store the mean square error of each cross validation
            mean_mse_list = []

            # Train and evaluate the model using K-fold cross validation
            for train_index, test_index in k_folds.split(X):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]

                # Create and train model
                model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=2024)
                model.fit(X_train, y_train)

                # Predict on the test set
                y_pred = model.predict(X_test)

                # Calculate mean square error
                mse = mean_squared_error(y_test, y_pred)
                mean_mse_list.append(mse)

            # Calculate the average of cross validation mean square error
            mean_mse = np.mean(mean_mse_list)

            # Update the best parameters and minimum mean square error
            if mean_mse < best_mean_mse:
                best_mean_mse = mean_mse
                best_alpha = alpha
                best_l1_ratio = l1_ratio

    #           print(f"alpha = {round(alpha, 2)}, l1_ratio = {round(l1_ratio, 3)}, mean MSE = {round(mean_mse, 3)}")

    # Output the best parameters
    print(f"Best alpha: {best_alpha}, Best l1_ratio: {best_l1_ratio}, Best mean MSE: {best_mean_mse}")
    t2 = time.time()
    print(f'10-fold cross validation is complete,takes {str(t2-t1)}')

    # Get 500 models through bootsrtap and get coefficients of elastic net regression
    print('Now is bootstrap')
    t3 = time.time()
    coef = bootstrap(bootstrap_data,500)
    coef = pd.DataFrame(coef,columns=list(data))
    t4 = time.time()
    print(f'bootstrap is complete,takes {str(t4-t3)}')

    # Delete features selected in less than half of the models (250)
    zero_percent = (coef == 0).sum(axis = 0) / 500
    columns_to_drop = zero_percent[zero_percent > 0.5].index

    # Output coefficients of important features in 500 models after screening
    coef_drop = coef.drop(columns_to_drop,axis = 1)

    # Determine the data to be used in the final model
    model_data = bootstrap_data.drop(columns_to_drop ,axis=1)
    print(f'The shape of model_data of {tissue} is: {model_data.shape}')
    model_data.to_csv(f'train_data4/bootstrap_{tissue}_{corr}.csv') 
```

## 3. train model

### 3.1 without disease features

```python
# Define parameter search range
alpha_range = [0.1,1,10,50]
l1_ratio_range = np.arange(0, 1, 0.05)

corr = 0.2
tissue_list = ['Blood','Brain','Lung','Skin']
for tissue in tissue_list:
    model_data = pd.read_csv(f'train_data4/bootstrap_{tissue}_{corr}.csv',index_col=0)

    loo = LeaveOneOut()
    # Store predicted age
    age_pred_list = []
    # Split the dataset using the leave-one-out method
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

        # K-fold cross validation to select the optimal parameters
        for alpha in alpha_range:
            for l1_ratio in l1_ratio_range:
                # Define K-fold cross validation
                k_folds = KFold(n_splits=5, shuffle=True, random_state=2024)

                # Store the mean square error of each cross validation
                mean_mse_list = []
                # Train and evaluate the model using K-fold cross validation
                for train_index, test_index in k_folds.split(train_x):
                    X_train, X_test = train_x.iloc[train_index], train_x.iloc[test_index]
                    y_train, y_test = train_y.iloc[train_index], train_y.iloc[test_index]

                    # Create and train model
                    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=2024)
                    model.fit(X_train, y_train)

                    # Predict on the test set
                    y_pred = model.predict(X_test)

                    # Calculate mean square error
                    mse = mean_squared_error(y_test, y_pred)
                    mean_mse_list.append(mse)

                # Calculate the average of cross validation mean square error
                mean_mse = np.mean(mean_mse_list)

                # Update the best parameters and minimum mean square error
                if mean_mse < best_mean_mse:
                    best_mean_mse = mean_mse
                    best_alpha = alpha
                    best_l1_ratio = l1_ratio

#                     print(f"alpha = {round(alpha, 2)}, l1_ratio = {round(l1_ratio, 3)}, mean MSE = {round(mean_mse, 3)}")

        # Output the best parameters
        print(f"Best alpha: {best_alpha}, Best l1_ratio: {best_l1_ratio}, Best mean MSE: {best_mean_mse}")

        # Leave-one-out method to predict age
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
    # Rearrange samples according to meta data
    model_data = model_data.loc[meta['SRA']]
    # Add columns based on disease information
    for i in set(meta['disease state']):
        model_data[i] = (meta['disease state']==i).astype(int).tolist()
    # Adjust the age column to be the last column
    model_data = model_data[[col for col in model_data.columns if col != 'age'] + ['age']]
    # save data
    model_data.to_csv(f'train_data4/bootstrap_{tissue}_{corr}_add.csv') 
```

```
# train model
# Define parameter search range
alpha_range = [0.1,1,10,50]
l1_ratio_range = np.arange(0, 1, 0.05)

corr = 0.2
tissue_list = ['Blood','Brain','Lung','Skin']
for tissue in tissue_list:
    model_data = pd.read_csv(f'train_data4/bootstrap_{tissue}_{corr}_add.csv',index_col=0)

    loo = LeaveOneOut()
    # Store predicted age
    age_pred_list = []
    # Split the dataset using the leave-one-out method
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

        # K-fold cross validation to select the optimal parameters
        for alpha in alpha_range:
            for l1_ratio in l1_ratio_range:
                # Define K-fold cross validation
                k_folds = KFold(n_splits=5, shuffle=True, random_state=2024)

                # Store the mean square error of each cross validation
                mean_mse_list = []
                # Train and evaluate the model using K-fold cross validation
                for train_index, test_index in k_folds.split(train_x):
                    X_train, X_test = train_x.iloc[train_index], train_x.iloc[test_index]
                    y_train, y_test = train_y.iloc[train_index], train_y.iloc[test_index]

                    # Create and train model
                    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=2024)
                    model.fit(X_train, y_train)

                    # Predict on the test set
                    y_pred = model.predict(X_test)

                    # Calculate mean square error
                    mse = mean_squared_error(y_test, y_pred)
                    mean_mse_list.append(mse)

                # Calculate the average of cross validation mean square error
                mean_mse = np.mean(mean_mse_list)

                # Update the best parameters and minimum mean square error
                if mean_mse < best_mean_mse:
                    best_mean_mse = mean_mse
                    best_alpha = alpha
                    best_l1_ratio = l1_ratio

#                     print(f"alpha = {round(alpha, 2)}, l1_ratio = {round(l1_ratio, 3)}, mean MSE = {round(mean_mse, 3)}")

        # Output the best parameters
        print(f"Best alpha: {best_alpha}, Best l1_ratio: {best_l1_ratio}, Best mean MSE: {best_mean_mse}")

        # Leave-one-out method to predict age
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

