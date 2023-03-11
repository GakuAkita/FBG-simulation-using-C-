# FBG Simulation

### モード結合理論を用いたFBGの反射スペクトルの計算に用いたスクリプトです。<br>
pythonは基本的にfor文が非常に遅く、numpyを使ってfor文を使わずにかければそれで良かったですが、計算過程が複雑でできませんでした。<br>
そこで、計算処理が早いC++を用いて行列の計算は行い、計算結果をcsvで出力しグラフ化はpythonで行うというような方法でシミュレーションしました。C++はfor文やそれだけに限らず処理が速いそうです。<br>
<br>
C++で高速に計算を行う方法を探す中で、Eigenというライブラリの存在を知りました。具体的に使うにはフォルダをダウンロードする必要があるが、それも含めてセットアップの手順は[この動画](https://www.youtube.com/watch?v=XmtNr1TuO-E)が参考になりました<br>
Source.cppがVisual Studioを用いて私が書いたコードであるが、こちらを参考にしながら0から手を動かしながら自分なりに書いていくことをおすすめします。
何か質問があれば気軽に連絡してください。
