////////////////////////////////////////////////////////////////
//                                                            //
// Fernanda Cavalcante Nascimento - 11390827                  //
//                                                            //
////////////////////////////////////////////////////////////////

import java.io.*;
import java.util.*;

class Vector {

	private double x;	// coordenada horizontal
	private double y;	// coordenada vertical

	public Vector(double x, double y){

		this.x = x;
		this.y = y;
	}

	public double getX() {

		return x;
	}

	public double getY() {
		
		return y;
	}

	// metodo que converte as coordenadas de um vetor do espaco referente a area de desenho
	// (veja detalhamento nos comentarios da funcao draw_line) para o espaco referente a imagem
	// propriamente dita.

	Vector convert(double width, double height){

		return new Vector((this.x + 1) * width / 2.0, (1 - this.y) * height / 2.0);
	}
}	

class Shape {

	public static final int MAX_VERTICES = 100;	// numero maximo de vertices que um shape pode possuir

	private int free;
	private Vector [] vertices;

	public Shape(int n){

		free = 0;
		vertices = new Vector[n];
	}

	public void add(Vector v){

		if(free < vertices.length){

			vertices[free++] = v;
		}
	}

	public int nVertices(){

		return vertices.length
;
	}

	public Vector get(int i){

		return i < free ? vertices[i] : null;
		
	}

}

class Image {

	public static final int MAX_VALUE = 255;	// valor maximo que um pixel da image pode assumir

	private int width;	// largura da imagem (numero de colunas da imagem)
	private int height;	// altura da imagem (numero de linhas da imagem)
	int [][] matrix;	// matriz que armazena o conteudo da imagem

	public Image(int width, int height, int bg_color){

		this.width = width;
		this.height = height;

		matrix = new int[this.height][this.width];

		for(int i = 0; i < this.height; i++) {

			for(int j = 0; j < this.width; j++) {

				matrix[i][j] = bg_color;
			}
		}
	}

	public int getWidth() { 
		
		return width; 
	}

	public int getHeight() { 

		return height;
	}

	// metodo que salva uma imagem em um arquivo, no formato PGM. Programas/utilitarios que trabalham com
	// imagens bitmap (como o Gimp e visualizadores de imagens normalmente instalados em ambiente Linux) em
	// geral reconhecem este formato. Se tiver dificuldades para visualizar os arquivos de imagem salvos por 
	// esta funcao, um visualizador online pode ser encontrado neste link: http://paulcuth.me.uk/netpbm-viewer/.

	public void save(String file_name){

		try{
			PrintWriter out = new PrintWriter(file_name);

			out.println("P2");
			out.println(width + " " + height);
			out.println(MAX_VALUE);
	
			for(int lin = 0; lin < height; lin++){

				for(int col = 0; col < width; col++){

					out.print((col == 0 ? "" : " ") + matrix[lin][col]); 				
				}
		
				out.println();
			}

			out.close();
		}
		catch(IOException e){
		
			e.printStackTrace();		
		}
	}

	// metodo que pinta um pixel da imagem com a "cor" (tonalidade de cinza) especificada.

	public void set_pixel(double x, double y, int color){

		int col = (int) Math.round(x);
		int lin = (int) Math.round(y);

		if(lin < 0 || col < 0 || col >= width || lin >= height) return;

		matrix[lin][col] = color;
	}

	// metodo que desenha uma linha conectando os pontos representados pelos vetores v1 e v2, 
	// com a cor especificada. A area de desenho compreende o retangulo formado pelos seguintes 
	// cantos, independente da dimensao real da imagem em pixels:
	//
	//   (-1,  1): canto superior esquerdo
	//   ( 1,  1): canto superior direito
	//   ( 1, -1): canto inferior direito
	//   (-1, -1): canto inferior esquerdo
	//
	// Logo, espera-se que as coordenadas dos vetores estejam dentro destes limites

	public void draw_line(Vector v1, Vector v2, int color){

		Vector p1 = v1.convert(width, height);
		Vector p2 = v2.convert(width, height);
		Vector a, b;

		double deltaX = Math.abs(p1.getX() - p2.getX());
		double deltaY = Math.abs(p1.getY() - p2.getY());
		double x, y;
	
		if(deltaX >= deltaY){

			a = p1.getX() < p2.getX() ? p1 : p2;
			b = p1.getX() < p2.getX() ? p2 : p1;

			for(x = a.getX(); x <= b.getX(); x++){

				y = ((x - a.getX()) / deltaX ) * (b.getY() - a.getY()) + a.getY();
				set_pixel(x, y, color);
			}  					
		}
		else{

			a = p1.getY() < p2.getY() ? p1 : p2;
			b = p1.getY() < p2.getY() ? p2 : p1;

			for(y = a.getY(); y <= b.getY(); y++){

				x = ((y - a.getY()) / deltaY ) * (b.getX() - a.getX()) + a.getX();
				set_pixel(x, y, color);
			}
		}
	}
}

class Matrix {

	public static final double SMALL = 0.000001;	// constante usada na comparação entre valores double

	private int lin, col;
	private double [][] m;

	// metodo estatico que cria uma matriz identidade de tamanho n x n. 

	public static Matrix identity(int n){

		Matrix mat = new Matrix(n, n);

		for(int i = 0; i < mat.lin; i++) mat.m[i][i] = 1.0;

		return mat;
	}

	// constroi uma matriz com as entradas iguais a zero

	public Matrix(int lin, int col){

		this.lin = lin;
		this.col = col;
		this.m = new double[this.lin][this.col];
	}

	// constroi uma matriz e inicia as entradas a partir do vetor values.

	public Matrix(int lin, int col, double [] values){

		int k = 0;
		
		this.lin = lin;
		this.col = col;
		this.m = new double[this.lin][this.col];		

		for(int i = 0; i < this.lin; i++){

			for(int j = 0; j < this.col; j++){

				this.m[i][j] = values[k++];
			}
		}
	}

	// metodo que devolve uma copia da matriz.

	Matrix copy(){
	
		Matrix copy = new Matrix(this.lin, this.col);

		for(int i = 0; i < this.lin; i++){

			for(int j = 0; j < this.col; j++){

				copy.m[i][j] = this.m[i][j];
			}
		}

		return copy;
	}

	// metodo que imprime a matriz.

	void print(){

		for(int i = 0; i < this.lin; i++){

			for(int j = 0; j < this.col; j++){
	
				System.out.printf("%7.2f ", this.m[i][j]);
			}

			System.out.println();
		}
	}

	// devolve o valor da entrada na linha i e coluna j.

	public double get(int i, int j){

		return this.m[i][j];
	}

	// define um novo valor para a entrada na coluna i e linha j.

	public void set(int i, int j, double value){

		this.m[i][j] = value;
	}

	// metodo que calcula o determinante de uma matriz 2x2.

	public double det2x2(){

		if(this.lin == 2 && this.col == 2){

			return this.m[0][0] * this.m[1][1] - this.m[0][1] * this.m[1][0];
		}

		throw new IllegalStateException("Invalid matrix size!");
	}

	// metodo que calcula o determinante de uma matriz 3x3.

	double det3x3(){

		if(this.lin == 3 && this.col == 3){

			return this.m[0][0] * this.cofactor3x3(0, 0) + this.m[0][1] * this.cofactor3x3(0, 1) + this.m[0][2] * this.cofactor3x3(0, 2);
		}

		throw new IllegalStateException("Invalid matrix size!");
	}

	// metodo que calcula o cofator para os indices (i, j) de uma matriz 3x3.

	public double cofactor3x3(int i, int j){

		if(this.lin == 3 && this.col == 3){

			int k = 0;
			double [] values = new double[4];
	
			for(int lin = 0; lin < this.lin; lin++){

				for(int col = 0; col < this.col; col++){

					if(lin != i && col != j) values[k++] = this.m[lin][col];
				}
			}

			Matrix tmp = new Matrix(2, 2, values);
	
			double cof = ((i + j) % 2 == 0 ? 1 : -1) * tmp.det2x2();

			return cof;

		}

		throw new IllegalStateException("Invalid matrix size!");
	}

	// funcao que calcula a inversa de uma matriz 3x3.

	public Matrix invert3x3(){

		double det = this.det3x3();

		if(Math.abs(det) < SMALL){

			throw new IllegalStateException("Singular matrix!\n");
		}

		Matrix inverse = new Matrix(3, 3);

		for(int i = 0; i < inverse.lin; i++){

			for(int j = 0; j < inverse.col; j++){

				inverse.m[j][i] = this.cofactor3x3(i, j) / det;
			}
		}

		return inverse;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                                                 //
	// Funcoes relacionadas a classe Matrix que precisam ser implementadas para o programa funcionar da forma esperada //
	//                                                                                                                 //
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	public Matrix multiply(Matrix m){

        int M1L=a->lin, M1C=a->col, M2L=b->lin, M2C=b->col;
        for(linha=0; linha<M1L; linha++){
            for(coluna=0; coluna<M2C; coluna++){
                somaprod=0;
                for(i=0; i<M1L; i++){
                    somaprod+=a->m[linha][i]*b->m[i][coluna];
                }
                resultante->m[linha][coluna]=somaprod;
            }
        }
        return resultante;
            
            return null;
    }

	public Vector transform(Vector v){

		// TODO: implementar!
		
		return v;
	}

	public static Matrix get_rotation_matrix(double theta){

		// TODO: implementar!
		
		return Matrix.identity(3);
	}

	public static Matrix get_scale_matrix(double k){

		// TODO: implementar!
		
		return Matrix.identity(3);
	}

	public static Matrix get_translation_matrix(Vector v){

		// TODO: implementar!
		
		return Matrix.identity(3);
	}

	public static Matrix get_transformation_matrix(Vector e1, Vector e2, Vector t){

		// TODO: implementar!
		
		return Matrix.identity(3);
	}

	public static Matrix get_observer_matrix(Vector position, Vector direction){

		// TODO: implementar!
		
		return Matrix.identity(3);
	}
}

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Classe principal. Le um arquivo de entrada com a seguinte estrutura:   //
//                                                                        //
//	IMAGE_WIDTH IMAGE_HEIGHT BG_COLOR                                 //
//	OBSERVER_X, OBSERVER_Y, DIRECTION_X, DIRECTION_Y                  //
//	N_SHAPES                                                          //
//	N_VERTICES_SHAPE0 X_0 Y_0 X_1 Y_1 ... X_(N-1) Y_(N-1)             //
//	N_VERTICES_SHAPE1 X_0 Y_0 X_1 Y_1 ... X_(N-1) Y_(N-1)             //
//	N_VERTICES_SHAPE2 X_0 Y_0 X_1 Y_1 ... X_(N-1) Y_(N-1)             //
//	...                                                               //
//	<DRAW_COMMAND_0>                                                  //
//	<DRAW_COMMAND_1>                                                  //
//	<DRAW_COMMAND_2>                                                  //
//	...                                                               //
//	END                                                               //
//                                                                        //
// Sendo que cada linha referente a um comando de desenho pode ser:       //
//                                                                        //
//	DRAW_SHAPE SHAPE_ID COLOR ROTATION SCALE T_X T_Y                  //
//                                                                        //
// OU                                                                     //
//                                                                        //
//	DRAW_SHAPE_BASE SHAPE_ID COLOR E1_X E1_Y E2_X E2_Y T_X T_Y        //
//                                                                        //
// E gera uma imagem a partir das configurações e comandos especificados. //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

public class EP2_esqueleto {

	public static String DRAW_SHAPE = "DRAW_SHAPE";			// constante para o comando DRAW_SHAPE usada no arquivo de entrada 
	public static String DRAW_SHAPE_BASE = "DRAW_SHAPE_BASE";	// constante para o comando DRAW_SHAPE_BASE usada no arquivo de entrada
	public static String END = "END";				// indica fim da listagem de comandos no arquivo de entrada

	public static void main(String [] args) throws IOException {

		////////////////
		//            //
		// Variaveis: //
		//            //
		////////////////

		String command;	
		String input_file_name;		// nome do arquivo de entrada (com as definicoes do cena/desenho)
		String output_file_name;	// nome do arquivo de saida (imagem gerada a partir das definicoes do arquivo de entrada)
	
		int width;			// largura da imagem a ser gerada
		int height;			// altura da imagem a ser gerada
		int background_color;		// cor de fundo da imagem a ser gerada
	
		int n_shapes;			// quantidade de shapes definidos no arquivo de entrada
	
		Vector observer;		// vetor que representa a posicao do observador na cena a ser desenhada
		Vector direction;		// vetor que indica a direcao para a qual o observador olha na cena
		Vector v;			// variavel do tipo vetor de uso geral
	
		Shape [] shapes;			// vetor de shapes. Usado para armazenar os shapes definidos no arquivo de entrada

		Scanner in;
		//FILE * in;			// FILE handler para o arquivo de entrada
		Image img;			// imagem na qual as operacoes de desenho serao realizada

		///////////////////////////////////////////
		//                                       //
		// Programa principal propriamente dito: //
		//                                       //
		///////////////////////////////////////////
	
		// verificacao dos parametros da linha de comando:
	
		if(args.length != 2){

			System.out.println("Usage: " + EP2_esqueleto.class.getName() + " <input_file_name> <output_file_name>");
			return;
		}
		
		input_file_name = args[0]; 
		output_file_name = args[1];
	
		// abertura do arquivo de entrada, e leitura dos parametros fixos (parametros da imagem e do observador, quantidade de shapes):

		in = new Scanner(new FileInputStream(input_file_name));

		width = in.nextInt();
		height = in.nextInt();
		background_color = in.nextInt();

		observer = new Vector(in.nextDouble(), in.nextDouble());
		direction = new Vector(in.nextDouble(), in.nextDouble());

		n_shapes = in.nextInt();
	
		img = new Image(width, height, background_color);	// criacao da imagem
		shapes = new Shape[n_shapes];				// alocacao do vetor de shapes com o tamanho adequado
	
		// leitura dos shapes definidos no arquivo de entrada:

		for(int i = 0; i < n_shapes; i++){
		
			shapes[i] = new Shape(in.nextInt());
	
			for(int j = 0; j < shapes[i].nVertices(); j++) {

				shapes[i].add(new Vector(in.nextDouble(), in.nextDouble()));
			}
		}

		// leitura dos comandos de desenho, até que uma linha com o comando "END" seja encontrada:

		while( !(command = in.next()).equals(END) ) {

			// a cada iteracao um shape sera desenhado atraves do comando DRAW_SHAPE ou DRAW_SHAPE_BASE.
			//		
			// Para primeiro comando (DRAW_SHAPE), especifica-se o id do shape a ser desenhdo, sua cor (tonalidade 
			// de cinza na realidade), assim como rotacao, fator de escala e vetor de translacao que devem ser 
			// aplicados ao shape antes de desenha-lo.
			//
			// Já para o segundo caso (DRAW_SHAPE_BASE), ao inves de especificar os fatores de rotacao e escala, 
			// especifica-se os vetores da base que define a transformacao matricial a ser aplicada no shape, além
			// de também se espeficiar o vetor de translacao que tambem deve ser aplicado. Apesar de ser menos intuitivo
			// este segundo comando permite especificar transformacoes com maior flexibilidade (escalas nao uniformes, e
			// cizalhamentos, por exemplo). 

			if(command.equals(DRAW_SHAPE) || command.equals(DRAW_SHAPE_BASE)){

				int shape_id = in.nextInt();
				int color = in.nextInt();;
			
				if(command.equals(DRAW_SHAPE)){

					double rotation = in.nextDouble();
					double scale = in.nextDouble();
					Vector t = new Vector(in.nextDouble(), in.nextDouble());

					// TODO: fala algo por aqui...
				}
			
				if(command.equals(DRAW_SHAPE_BASE)){

					Vector e1 = new Vector(in.nextDouble(), in.nextDouble());
					Vector e2 = new Vector(in.nextDouble(), in.nextDouble());
					Vector t = new Vector(in.nextDouble(), in.nextDouble());
			
					// TODO: fala algo por aqui...
				}

				Shape s = shapes[shape_id];

				for(int i = 0; i < s.nVertices() - 1; i++){

					// TODO: fala algo por aqui...

					Vector v1 = s.get(i);
					Vector v2 = s.get(i + 1);
					img.draw_line(v1, v2, color);
				}
			}
			else {
				System.out.println("Unknown command: " + command);
			}
		}

		img.save(output_file_name);	// salva imagem no arquivo de saida	
	}
}
