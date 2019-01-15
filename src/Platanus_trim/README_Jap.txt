Platanus_trim v1.0.3 マニュアル



1. Platanus_trim
1.1 概要
  以下の3つのトリムを上から順に行うプログラムです。
  1 read中のアダプター配列を除去
  2 pair同士の突き抜け(forwardとreverseで共通する部分をシーケンスしている箇所)を特定し除去
  3 qualityが低い領域と、上記1,2のトリムの結果極端に短くなったリード(11bp未満)の除去
  アダプター配列とは完全一致している必要はなく、ある程度のエラーを許容します。


1.2 使用方法
  $ ./platanus_trim [option] PAIR_FWD1 PAIR_REV1 (PAIR_FWD2 PAIR_REV2 ...)

  PAIR_FWD1, PAIR_REV1はペアのreadファイルです。
  1つのファイルにforwardとreverseが両方含まれている場合は分割してください。


1.3 オプション
  -i STR
    STRにトリムするリードのファイル一覧が記載されたファイルを指定することで
    逐次ファイルを指定せずにトリムを行うことができます。

  -q INT
    qualityトリムの閾値です。指定した値より大きいquality値が11base以上
    連続して続いたとき、その領域の先頭より前のreadをトリムします。

  -t INT
    本プログラムで使用するスレッド数を指定します。

  -1 STR
    2つあるアダプタ配列の1つをSTRで指定します。

  -2 STR
    2つあるアダプタ配列のもう1つをSTRで指定します。

    -1と-2は基本的には順序関係はありません
    どちらにどのアダプタ配列を書いても同じ動作をします。
    しかし、reverse complement考慮してアダプタ同士に共通の11merが存在した場合
    -1で指定した側のアダプタが上手く取れない不具合があります。
    この場合は、トリムされた配列を入力とし-1と-2のオプションを入れ替えて
    再度実行してください。


1.4 アウトプット
  トリムされたreadファイルが、入力ファイルが存在する同じディレクトリに
  (入力ファイル名).trimmed
  という名前で作成されます。


1.5 ログ
  アダプター配列の除去と突き抜けの除去、クオリティカットの2つのトリムの結果として
  以下のようなログが標準エラーに出力されます。

  ex.
    NUM_OF_TRIMMED_READ(FORWARD) = 32966
    NUM_OF_TRIMMED_BASE(FORWARD) = 2179910
    NUM_OF_TRIMMED_READ(REVERSE) = 53178
    NUM_OF_TRIMMED_BASE(REVERSE) = 1850338
    NUM_OF_TRIMMED_PAIR(OR) = 78649
    NUM_OF_TRIMMED_PAIR(AND) = 7495

  これらの数値は上から順に
    実行時に指定したpairファイルのうち、先に指定したファイル中にあるreadのうちトリムされた本数
    実行時に指定したpairファイルのうち、先に指定したファイル中にあるreadのうちトリムされた塩基数
    実行時に指定したpairファイルのうち、後に指定したファイル中にあるreadのうちトリムされた本数
    実行時に指定したpairファイルのうち、後に指定したファイル中にあるreadのうちトリムされた塩基数
    readをペア単位で見たとき、少なくとも片方のreadがトリムされたペア数
    readをペア単位で見たとき、両方のreadがトリムされたペア数
  を表しています。








2. Platanus_internal_trim
2.1 概要
  mate-pair作成時に挿入されることがあるinternal adaptorを除去するプログラムです。
  原則perfect matchですが、readの末端にアダプターの一部分だけ含まれている場合には
  18mer以上のperfect matchがあると除去します。
  このプログラムではplatanus_trimの機能も有しており以下の挙動を示します。
  1 internal adaptorのトリムを行う。
  2 internal adaptorがなかったreadに対してplatanus_trimと同じ動作をさせる。


1.2 使用方法
  $ ./platanus_internal_trim [option] PAIR_FWD1 PAIR_REV1 (PAIR_FWD2 PAIR_REV2 ...)

  PAIR_FWD1, PAIR_REV1はペアのreadファイルです。
  1つのファイルにforwardとreverseが両方含まれている場合は分割してください。


1.3 オプション
  -a INT
    internal adaptorの配列を指定します。
    INTの数値と対応するアダプターについてはplatanus_internal_trimを引数なしで実行してください。

  -b STR
    internal adaptorを手動で指定します。
    STRにinternal adaptorの配列を指定します。
    なお、-aと同時に使用した場合-aが優先されます。

  -r
    internal
    adaptorのカットのみを行い、platanus_trimに相当する動作を行わないようにします。

  またplatanus_trimにある全てのオプションが使用でき、その動作についてはplatanus_trimと全く同じです。


1.4 アウトプット
  トリムされたreadファイルが、入力ファイルが存在する同じディレクトリに
  (入力ファイル名).int_trimmed
  という名前で作成されます。


1.5 ログ
  platanus_trimのログに加えて、internal adaptorによってトリムされたログが追加されます。
  詳細は1.5を参照のこと。









3. Tips
3.1 platanus_trim(platanus_internal_trim)でクオリティカットしたくない
  -qオプションに0を指定してください。


3.2 platanus_trim(platanus_internal_trim)でクオリティカットだけしたい
  -1と-2オプションに空文字を指定してください。
  ex. $ ./platanus_trim -1 "" -2 ""


3.4 platanus_trim(platanus_internal_trim)で-q 0にも関わらず
ログのクオリティでカットされたread数が1以上の時がある
  クオリティカットは、アダプターや突き抜けのトリムを行った後のreadが11bpより小さいものを
  自動的にトリムする機能が付いています。
  この機能によってトリムされたreadのカウントが行われています。







Author
Kota Toshimoto
