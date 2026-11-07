# SEA Open

リガンド集合の化学的類似性から標的を関連付けるSimilarity Ensemble Approach (SEA) のOSS再実装です。Pixi、RDKit、FastAPIを使用します。日本語Web UI、APIキー認証、永続ジョブ、独立workerを備えています。

原著: Keiser et al., *Relating protein pharmacology by ligand chemistry*, Nature Biotechnology 25, 197–206 (2007), https://doi.org/10.1038/nbt1284

## 起動

Linux x86_64とPixiを使用します。

```bash
pixi install --locked
pixi run demo
pixi run sea key alice
```

最後のコマンドが表示するAPIキーを保存してください。キーはDBにはSHA256として保存されます。別々のターミナルで実行します。

```bash
pixi run worker
pixi run serve
```

http://127.0.0.1:8000 を開き、APIキーで接続します。API文書は `/docs`、機械可読仕様は `/openapi.json` です。デモは合成のアルコール・カルボン酸・芳香族集合であり、生物学的標的の予測に使用できません。

`SEA_HOME`（既定 `var/`）に参照DBとサービスDBを保存します。`pixi run sea datasets` でdataset IDを取得します。サンドボックスで `/tmp` に書けない場合は、作業ディレクトリ内に一時ディレクトリを作り `TMPDIR="$PWD/.cache/tmp" pixi …` を使用してください。

## ChEMBL参照DB

[公式配布](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/)から、使用するリリースと対応するSHA256を確認します。版を固定し、数GBの取得と展開領域を確保してください。自動で最新版へ更新しません。

```bash
pixi run sea download --release RELEASE_NUMBER --sha256 OFFICIAL_SHA256
pixi run sea build --sqlite var/downloads/chembl_RELEASE_NUMBER.db --release ChEMBL_RELEASE_NUMBER
```

`RELEASE_NUMBER`、`OFFICIAL_SHA256` は実際の値に置換します。既存の公式SQLiteを `--sqlite` に直接指定することもできます。ファイル内のversionテーブルと指定版を照合します。

既定条件は `config/chembl.json`：ヒト、単一タンパク質、confidence 9、binding assay、Ki/Kd/IC50、厳密な等号測定、有効なpChEMBL、品質問題・potential duplicate除外。親分子構造を正規化し、標的・構造・測定指標ごとの中央値がpChEMBL ≥ 6のものを活性リガンドとします。10分子以上の標的を採用します。指標間の平均は計算しません。

標的、リガンド、活性測定の由来、原データSHA256、フィルター、除外件数、RDKit版、校正情報を不変のJSONに保存します。設定・データ・校正が変わると新しいdataset IDになります。ソースデータと校正を変更しない再構築では同じIDになります。

```bash
pixi run sea calibrate DATASET_ID --sizes 1,10,30,100,300,1000 --repeats 128
pixi run sea evaluate DATASET_ID --output var/evaluation.json
```

大きい集合の校正は計算時間が増えます。初回構築は最大100分子で校正し、それ以上はMonte Carloで評価します。評価コマンドは化合物を全標的横断のMurcko scaffoldで分割し、学習側だけで背景モデルを再推定します。既知標的へのRecall@1/5/10・MRRを最大Tanimoto法と比較します。未測定ペアを陰性とは見なしません。既定で最大100クエリを評価します。

## Agentからの利用

環境変数 `SEA_API_KEY` にキー、`SEA_DATASET_ID` にdataset IDを設定します。

```bash
curl -H "X-API-Key: $SEA_API_KEY" http://127.0.0.1:8000/api/v1/datasets
```

[examples/agent.py](examples/agent.py) は要求の送信、状態確認、結果取得までを実行します。

```bash
pixi run python examples/agent.py
```

`POST /api/v1/predictions` の `queries` は集合の配列です。`[["CCO"],["CCCO"]]` は分子別検索、`[["CCO","CCCO"]]` は1集合の検索です。`POST /api/v1/comparisons` は `left` と `right` を受け取ります。どちらも `dataset_id` が必須で、HTTP 202と `job_id` を返します。

- `GET /api/v1/jobs/{id}`：状態と進捗
- `GET /api/v1/jobs/{id}/result`：JSON、`?format=csv` でCSV
- `POST /api/v1/jobs/{id}/cancel`：取消
- `GET /api/v1/jobs`：自分の直近100件
- `GET /api/v1/targets?dataset_id=...&q=...`：標的検索

既定上限は要求あたり100分子、集合比較は片側1,000分子、body 10 MB。`SEA_MAX_QUERY_MOLECULES`、`SEA_MAX_SET_MOLECULES`、`SEA_MAX_UPLOAD_BYTES`で変更できます。入力不正はHTTP 422、認証失敗401、他ユーザーのジョブは404、未完了の結果取得は409です。

## 研究室での共有

同一ホストでAPIと `pixi run worker` を実行します。workerは既定2プロセス、`pixi run sea worker --count 4` で変更可能です。supervisor再起動時は中断ジョブを再キューします。計算は決定的で、再実行による外部への書込みはありません。サービスDBとデータセットをバックアップしてください。ジョブは自動削除しません。

利用者ごとにキーを発行してください。ログインは8時間のHttpOnly・SameSite=Strict cookieを使用し、APIキーをブラウザストレージには保存しません。TLSを使用する共有環境では `SEA_SECURE_COOKIE=1` とし、[docs/Caddyfile](docs/Caddyfile)を参考にTLSリバースプロキシを配置します。APIキー発行・DB構築はサーバー管理者のCLI操作です。APIは127.0.0.1へバインドし、プロキシだけを研究室ネットワークへ公開します。

## 検証と制約

```bash
pixi run test
```

計算詳細は [docs/method.md](docs/method.md)、今回の動作検証は [docs/validation.md](docs/validation.md) を参照してください。原著のDaylightフィンガープリントとMDDRをRDKitとChEMBLへ置換するため、原著の数値一致は主張しません。E値は結合確率ではありません。Monte Carloの既定試行数999ではp値の最小値は0.001であり、多数の標的を検索する場合に有意性の分解能が不足することがあります。

ソースコードはMIT。合成デモはCC0。ChEMBL由来データはCC BY-SA 3.0で、ソフトウェアとは別のライセンスです。構築データの再配布時にはmanifestの出典・版情報も含めてください。

Monte Carlo試行数は `sea build --trials 99999` または `sea calibrate DATASET_ID --trials 99999` で変更できます。試行数を増やすとp値の分解能が上がり、計算時間も増えます。結果の `execution` にジョブ時間とworkerのピークRSS（KiB）を記録します。
