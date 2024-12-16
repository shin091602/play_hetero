# play_hetero

ヘテロクリニック接続の研究用のリポジトリです．

## commit message

- fix: バグ修正
- add: 新規機能追加
- updata: 機能修正
- remove: 削除

## Script

### script

- check_eigenvalue.m:
  Lyapunov orbit の family を計算し，その Monodoromy matrix の固有値を計算する
- detect_lyapunov.m:
  ターゲットとなるヤコビ定数の Lyapunov orbit を計算し，多様体を計算する
- tori_manifolds.m:
  トーラスの多様体を計算する

### function

- fun_manifolds_custom.m:　
  絶対値が１でない，面外方向の固有ベクトルの方向に多様体を計算する
- fun_manifolds_double.m:
  絶対値が１でない，２組の複素共役な固有値の固有ベクトルの方向に多様体を計算する
