# 化合物分割モジュール

Spressoで利用していた化合物分割モジュールを分離したものです。

## 環境構築

### Docker
.devcontainer/Dockerfile にコンテナ構築用のファイルが存在します。
Docker（特に VSCode のDocker extension）を利用されている方はこちらの利用が便利です。

コンテナ構築に一定の時間を要するためご注意ください（boostのインストールにそれなりに時間がかかります）。

### 通常インストール

#### 必要なライブラリ
本ツールには以下のライブラリが必要です。
- boost (version >= 1.36)
  - `sudo apt-get install libboost-all-dev`
- openbabel (version== 2.4.1)
  - `sudo apt-get install libopenbabel-dev`
  - **openbabel 3 では正常動作しません**。2.4.1 でのみ稼働します。

#### ビルド

`make` コマンドを実行すると実行可能ファイル `decompose` が作られます。

このビルドには、環境変数 `BOOST_INSTALL_PATH` および `OBABEL_INSTALL_PATH` を利用します。
もしこれら2つが `/usr` 以外に存在する場合（例えば手動でインストールした場合など）には、
環境変数を設定してください。

#### 使い方

コンフィグファイルを使う形式と、直接コマンドライン引数でファイルを指定する形式の2つに対応しています。
```
$ ./decompose -c testdata/decomposition_test.in
$ ./decompose -l testdata/G39.mol2 -f testdata/G39fragment.sdf -o testdata/G39annotated.sdf
```

なお、 `-h` オプションをつけると、以下のようなヘルプが表示されます。

```
$ ./decompose -h
Usage: ./decompose[options]
Options:
  -h [ --help ]                  show help
  -c [ --conf-file ] arg         configuration file
  -f [ --fragment ] arg          fragment file (.mol2 file)
  -l [ --ligand ] arg            ligand file
  -o [ --output ] arg            output (annotated) ligand file (.sdf file)
```

#### 入出力ファイルの説明

##### 入力ファイル

`-l [ --ligand ] LIGAND`
- 化合物分割を行う対象となる化合物ファイルです。
  - 1つのファイルに含まれる**複数の化合物を同時に処理**します。
  - **2次元構造か3次元構造を持つ化合物**を入力してください。
    - 座標がすべて0になっているデータでは正常に化合物分割が行えません。

##### 出力ファイル

` -f [ --fragment ] FRAGMENT`
- 化合物分割の結果作成されたフラグメントが保存されるファイルです。 `.sdf` 形式を推奨します。
  - 分子名はそのフラグメントのSMILES表記（OpenBabelのCanonical SMILES）になっています。
  - 各フラグメントデータは立体構造情報を持ちますが、これは入力されたファイルに含まれているデータに依存します。


` -o [ --output ] OUTPUT`
- 入力された化合物が、どのようなフラグメントからなるかが追記されたファイルです。 `.sdf` 形式を推奨します。
  - `-l` で入力された化合物情報はそのままに、 `fragment_info` というメタ情報が追加されています。
    - `fragment_info` は、その化合物が持つフラグメントのSMILES表記が、カンマ区切りで全て記載されています。

