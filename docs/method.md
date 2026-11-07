# 計算方法と原著との差分

1. RDKit Cleanup → FragmentParent → Uncharger → canonical isomeric SMILES。集合内の同一構造を除き、集合間の共通構造は残す。互変異性体は統合しない。
2. RDKit path fingerprint: 2048 bit, minPath=1, maxPath=7。他の設定はロックされたRDKit版の既定値。Tanimotoが閾値以上の全ペアのTanimotoを合計する（閾値を引き算しない）。
3. 校正用の背景分子は参照DB内のユニーク構造を一様にサンプリング。集合内は非復元抽出、2集合は独立に抽出する。集合間の重複は許す。
4. 各サイズ対で反復サンプルを生成し、前半で平均 `mu(n)=a*n`、標準偏差 `sigma(n)=b*n**c` を推定する。nは集合サイズの積。sigmaはlog空間で回帰する。
5. 閾値0.30–0.80（0.01刻み）でZ-scoreを計算し、Gumbel分布の位置・尺度を学習する。後半サンプルに対するKS統計量が最小の閾値を採用する。サイズ積100未満と分散ゼロは分布fitから除く。独立した最終評価はscaffold holdoutで行う。
6. held-out KS ≤ 0.10、両集合10分子以上、校正サイズ範囲内の場合にGumbelによるp値を使用。それ以外は実際のサイズ対でMonte Carloを行う。試行数999、`p=(exceed+1)/(trials+1)`。同じジョブ内でサイズ対の背景スコアを再利用する。集合サイズが背景プールを超える入力は拒否する。
7. `E=p*N`。Nは検索対象の全標的数（top_kではない）、2集合の比較では1。log pとlog Eは自然対数。Eは1を超え得る。順位はlog E、小さい順、同値ならraw scoreの大きい順。Monte Carloの分解能による同順位が生じる。

これは原著の統計的枠組みを使った実装上の選択であり、原著と同じ背景分布を保証しない。特に単一分子には連続極値分布を仮定せず、有限試行による評価を返す。閾値選択とGumbel適合の診断はdatasetのmodelに保存する。

原著の係数は専用テストにだけ使用する。原著Supplementary Methodsの数値例は丸め誤差があり、本文のE値との完全一致を期待しない。

- 原著: https://doi.org/10.1038/nbt1284
- Supplementary Methods: https://static-content.springer.com/esm/art%3A10.1038%2Fnbt1284/MediaObjects/41587_2007_BFnbt1284_MOESM64_ESM.pdf
- ChEMBL confidence: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions
- RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html

全標的ネットワークとMCP、一般公開向けアカウント管理、分散workerは初版に含まない。単一ホストのSQLiteキューは、研究室内の小規模共有向け。
