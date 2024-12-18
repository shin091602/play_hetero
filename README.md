# play_hetero

ヘテロクリニック接続の研究用のリポジトリです．

## commit message

- fix: バグ修正
- add: 新規機能追加
- updata: 機能修正
- remove: 削除

## script

- check_eigenvalue.m:
  Lyapunov orbit の family を計算し，その Monodoromy matrix の固有値を計算する
- detect_lyapunov.m:
  ターゲットとなるヤコビ定数の Lyapunov orbit を計算し，多様体を計算する
- tori_manifolds.m:
  トーラスの多様体を計算する
- hetero_tori_lyapunov.m:
  二次元トーラスの不安定多様体とリアプノフ軌道の安定多様体をそれぞれ伝播させ，接続を調べる

## function

- fun_manifolds_custom.m:　
  絶対値が１でない，面外方向の固有ベクトルの方向に多様体を計算する
- fun_manifolds_double.m:
  絶対値が１でない，２組の複素共役な固有値の固有ベクトルの方向に多様体を計算する

## Detail of Script

### hetero_tori_lyapunov.m

- まず，トーラスを計算し，そのヤコビ定数を調べる．その後同じヤコビ定数を持つリアプノフ軌道を計算し多様体を伝播させる
