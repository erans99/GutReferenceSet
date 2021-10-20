def simple_statistics(df):
    print()
    print('df shape:', df.shape[0])
    for col in df.columns:
        print()
        print(col)
        print('# missing values:', df[col].isna().sum())
        unique_values = set(df[col].dropna().unique())
        if 'Unknown' in unique_values:
            print('# unknown values:', (df[col] == 'Unknown').sum())
            unique_values.discard('Unknown')
        print('# unique values :', len(unique_values), unique_values if len(unique_values) <= 10 else '')
