# 動作確認記録

実行環境：Linux x86_64、Pixi 0.75.0、Python 3.12.16、RDKit 2026.09.2。依存関係はpixi.lockで固定。

最終結果：`pixi run test` は **10 passed in 2.67s**。`pixi install --locked --offline` も成功。

## 確認した内容

- SEA計算：構造正規化、塩除去、重複排除、類似度閾値の境界、集合比較の対称性。
- 原著Supplementary Methodsの数式：Z-scoreと極小p値、log値の精度。
- Monte Carlo：固定seedでの決定性、試行数によるp値分解能、集合サイズ上限。
- ChEMBLスキーマを再現したSQLite fixture：quality filter、活性中央値、標的採用、由来情報、版照合、再構築のID一致。
- API：認証、所有者分離、入力エラー、CSV出力、Cookieセッション、Origin制御、body上限。
- ジョブ：排他的claim、待機中・実行中の取消、supervisor再起動後の中断ジョブ回復。
- Chrome headless：実API・2プロセスworkerで、ログイン→2分子の標的検索→6行の結果表示。デスクトップ1280 px、モバイル390 pxのスクリーンショットを目視確認。

ブラウザでの確認用コマンド（開発時のみNodeとChromeが必要）：

```bash
pixi run python scripts/smoke.py
```

計測例（合成30分子・3集合、2分子別検索、Monte Carlo 999試行）：ジョブ計算0.064秒、workerピークRSS 88,016 KiB。デモ構築・サーバー／worker／Chromeの起動・検索・スクリーンショット取得を含むスモークチェックは約4.6秒。これは小規模デモの測定であり、ChEMBL全体での性能を示す値ではない。

## 未実施の検証

ChEMBL全件のダウンロード・参照DB構築、および実データによるscaffold holdoutの生物学的予測精度評価は、この実行では行っていない。構築・校正・評価用CLIは実装済み。公開前の科学的評価には実データで `sea build` と `sea evaluate` を実行する必要がある。合成デモの標的名は化学構造のグループ名であり、実際のタンパク質ではない。

この環境のサンドボックスではTestClientの内部イベントループ通信が停止するため、API／ブラウザの確認はローカル通信を許可した実行環境で行った。
