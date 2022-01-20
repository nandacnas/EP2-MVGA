////////////////////////////////////////////////////////////////
//                                                            //
// Fernanda Cavalcante Nascimento - 11390827                  //
//                                                            //
////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

////////////////
//            //
// Constantes //
//            //
////////////////


#define MAX_VERTICES 100			// numero maximo de vertices que cada Shape pode ter
#define MAX_VALUE 255				// valor maximo que um pixel da image pode assumir

#define DRAW_SHAPE "DRAW_SHAPE"			// constante para o comando DRAW_SHAPE usada no arquivo de entrada 
#define DRAW_SHAPE_BASE "DRAW_SHAPE_BASE"	// constante para o comando DRAW_SHAPE_BASE usada no arquivo de entrada
#define END "END"				// indica fim da listagem de comandos no arquivo de entrada

#define SMALL 0.000001				// constante para comparação entre valores double
#define EQUAL(a, b) (fabs(a - b) < SMALL)	// macro para comparar se dois valores double sao considerados iguais
#define ZERO(a) (fabs(a) < SMALL)		// macro para verificar se um valor double eh considerado igual a zero


//////////////////////////////////////////
//                                      //
// Estruturas consideradas neste codigo //
//                                      //
//////////////////////////////////////////


typedef struct {	// representa um vetor - ou ponto - bidimensional

	double x;	// coordenada horizontal
	double y;	// coordenada vertical

} Vector;		


typedef struct {	// representa uma forma geometrica que nada mais eh do que uma sequencia de vetores - ou pontos.

	int n;
	Vector vertices[MAX_VERTICES];

} Shape;


typedef struct {	// estrutura para representar uma imagem digital como uma matriz de valores inteiros

	int width;	// largura da imagem (numero de colunas da imagem)
	int height;	// altura da imagem (numero de linhas da imagem)
	int ** matrix;	// ponteiro para a matriz que armazena o conteudo da imagem

} Image;


typedef struct {	// estrutura para representar matrizes no sentido matematico
	
	double ** m;
	int lin, col;

} Matrix;


/////////////////////////////////////////
//                                     //
// Funcoes relacionadas ao tipo Vector //
//                                     //
/////////////////////////////////////////


// funcao auxiliar: cria um vetor a partir de coordenadas x e y do tipo double.

Vector create_vector(double x, double y){

	Vector v;

	v.x = x;
	v.y = y;

	return v;
}

// funcao que converte as coordenadas de um vetor do espaco referente a area de desenho
// (veja detalhamento nos comentarios da funcao draw_line) para o espaco referente a imagem
// propriamente dita.

Vector convert(Vector v, double width, double height){

	Vector u;

	u.x = (v.x + 1) * width / 2.0 ;
	u.y = (1 - v.y) * height / 2.0;

	return u;
}


////////////////////////////////////////
//                                    //
// Funcoes relacionadas ao tipo Image //
//                                    //
////////////////////////////////////////


// funcao que cria uma imagem vazia, com a cor de fundo especificada e a devolve.

Image * create_imagem(int w, int h, int bg_color){

	int lin, col;
	Image * image = (Image *) malloc(sizeof(Image));

	image->height = h;
	image->width = w;
	image->matrix = (int **) malloc(image->height * sizeof(int *));

	for(lin = 0; lin < image->height; lin++){

		image->matrix[lin] = (int *) malloc(image->width * sizeof(int));

		for(col = 0; col < image->width; col++) image->matrix[lin][col] = bg_color;
	}

	return image;
}

// funcao que libera os recursos de memoria associados a uma imagem.

void free_image(Image * image){

	int lin;

	if(image){

		for(lin = 0; lin < image->height; lin++) free(image->matrix[lin]);

		free(image->matrix);
		free(image);
	}
}

// funcao que salva uma imagem em um arquivo, no formato PGM. Programas/utilitarios que trabalham com
// imagens bitmap (como o Gimp e visualizadores de imagens normalmente instalados em ambiente Linux) em
// geral reconhecem este formato. Se tiver dificuldades para visualizar os arquivos de imagem salvos por 
// esta funcao, um visualizador online pode ser encontrado neste link: http://paulcuth.me.uk/netpbm-viewer/.

void save_image(Image * image, char * file_name){

	int lin, col;
	FILE * out = fopen(file_name, "w");

	if(out){

		fprintf(out, "P2\n%d %d\n%d\n", image->width, image->height, MAX_VALUE);

		
		for(lin = 0; lin < image->height; lin++){

			for(col = 0; col < image->width; col++){

				fprintf(out, col == 0 ? "%d" : " %d", image->matrix[lin][col]);  				
			}
			
			fprintf(out, "\n");
		}
	}

	fclose(out);
}

// funcao que pinta um pixel da imagem com a cor especificada.

void set_pixel(Image * image, double x, double y, int color){

	int col = (int) round(x);
	int lin = (int) round(y);

	if(lin < 0 || col < 0 || col >= image->width || lin >= image->height) return;

	image->matrix[lin][col] = color;
}

// funcao que desenha uma linha conectando os pontos representados pelos vetores v1 e v2, 
// com a cor especificada. A area de desenho compreende o retangulo formado pelos seguintes 
// cantos, independente da dimensao real da imagem em pixels:
//
//   (-1,  1): canto superior esquerdo
//   ( 1,  1): canto superior direito
//   ( 1, -1): canto inferior direito
//   (-1, -1): canto inferior esquerdo
//
// Logo, espera-se que as coordenadas dos vetores estejam dentro destes limites

void draw_line(Image * image, Vector v1, Vector v2, int color){

	Vector p1 = convert(v1, image->width, image->height);
	Vector p2 = convert(v2, image->width, image->height);
	Vector a, b;

	double deltaX = fabs(p1.x - p2.x);
	double deltaY = fabs(p1.y - p2.y);
	double x, y;
	
	if(deltaX >= deltaY){

		a = p1.x < p2.x ? p1 : p2;
		b = p1.x < p2.x ? p2 : p1;

		for(x = a.x; x <= b.x; x++){

			y = ((x - a.x) / deltaX ) * (b.y - a.y) + a.y;
			set_pixel(image, x, y, color);
		}  					
	}
	else{

		a = p1.y < p2.y ? p1 : p2;
		b = p1.y < p2.y ? p2 : p1;

		for(y = a.y; y <= b.y; y++){

			x = ((y - a.y) / deltaY ) * (b.x - a.x) + a.x;
			set_pixel(image, x, y, color);
		}
	}
}


/////////////////////////////////////////
//                                     //
// Funcoes relacionadas ao tipo Matrix //
//                                     //
/////////////////////////////////////////


// funcao que aloca uma matriz de n linhas por m colunas, e inicia as entradas com o valor zero.

Matrix * create_matrix(int n, int m){
	
	int i;

	Matrix * mat = (Matrix *) malloc (sizeof(Matrix)); 

	mat->lin = n;
	mat->col = m;
	mat->m = (double **) malloc(mat->lin * sizeof(double *));

	for(i = 0; i < mat->lin; i++){

		mat->m[i] = (double *) malloc(mat->col * sizeof(double));
		memset(mat->m[i], 0, mat->col * sizeof(double));
	}

	return mat;
}

// funcao que aloca uma matriz de n linhas por m colunas, e inicia as entradas a partir do vetor values.

Matrix * create_matrix_with_values(int n, int m, double * values){
	
	int i, j, k;

	Matrix * mat = create_matrix(n, m);

	k = 0;

	for(i = 0; i < mat->lin; i++){

		for(j = 0; j < mat->col; j++){

			mat->m[i][j] = values[k++];
		}
	}

	return mat;
}

// funcao que devolve uma copia da matriz recebida como parametro.

Matrix * copy_matrix(Matrix *  mat){
	
	int i, j;

	Matrix * copy = create_matrix(mat->lin, mat->col);

	for(i = 0; i < copy->lin; i++){

		for(j = 0; j < copy->col; j++){

			copy->m[i][j] = mat->m[i][j];
		}
	}

	return copy;
}

// funcao que cria uma matriz identidade com n linhas por n colunas. 

Matrix * create_identity(int n){

	int i;

	Matrix * mat = create_matrix(n, n);

	for(i = 0; i < mat->lin; i++) mat->m[i][i] = 1;

	return mat;
}

// funcao que imprime a matriz passada como parametro.

void print_matrix(Matrix * mat){

	int i, j;

	if(!mat) {

		printf("NULL MATRIX\n");
		return;
	}

	for(i = 0; i < mat->lin; i++){

		for(j = 0; j < mat->col; j++){
	
			printf("%7.2f ", mat->m[i][j]);
		}

		printf("\n");
	}
}

// funcao que libera os recursos de memoria associados a uma matriz.

void free_matrix(Matrix * mat){

	int lin;

	if(mat){

		for(lin = 0; lin < mat->lin; lin++) free(mat->m[lin]);

		free(mat->m);
		free(mat);
	}
}

// funcao que calcula o determinante de uma matriz 2x2.

double det2x2(Matrix * mat){

	if(mat->lin == 2 && mat->col == 2){

		return mat->m[0][0] * mat->m[1][1] - mat->m[0][1] * mat->m[1][0];
	}

	printf("Invalid matrix size!\n");
	return 0.0;
}

// funcao que calcula o cofator para os indices (i, j) de uma matriz 3x3.

double cofactor3x3(Matrix * mat, int i, int j){

	double cof;
	double values[4];
	int lin, col, k;	
	Matrix * tmp;

	if(mat->lin == 3 && mat->col == 3){

		k = 0;
	
		for(lin = 0; lin < mat->lin; lin++){

			for(col = 0; col < mat->col; col++){

				if(lin != i && col != j) values[k++] = mat->m[lin][col];
			}
		}

		tmp = create_matrix_with_values(2, 2, values);
	
		cof = ((i + j) % 2 == 0 ? 1 : -1) * det2x2(tmp);

		free(tmp);

		return cof;

	}

	printf("Invalid matrix size!\n");
	return 0.0;	
}

// funcao que calcula o determinante de uma matriz 3x3.

double det3x3(Matrix * mat){

	if(mat->lin == 3 && mat->col == 3){

		return mat->m[0][0] * cofactor3x3(mat, 0, 0) + mat->m[0][1] * cofactor3x3(mat, 0, 1) + mat->m[0][2] * cofactor3x3(mat, 0, 2);
	}

	printf("Invalid matrix size!\n");
	return 0.0;
}

// funcao que calcula a inversa de uma matriz 3x3.

Matrix * invert3x3(Matrix * mat){

	int i, j;
	double det = det3x3(mat);
	Matrix * inverse = create_matrix(3, 3);

	if(ZERO(det)){

		printf("singular matrix!\n");
		return NULL;
	}

	for(i = 0; i < inverse->lin; i++){

		for(j = 0; j < inverse->col; j++){

			inverse->m[j][i] = cofactor3x3(mat, i, j) / det;
		}
	}

	return inverse;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                //
// Funcoes relacionadas ao tipo Matrix que precisam ser implementadas para o programa funcionar da forma esperada //
//                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix * multiply(Matrix * a, Matrix * b){

	Matrix * resultante = (Matrix *) malloc (sizeof(Matrix)); 

	resultante->lin = a->lin;
	resultante->col = b->col;

	int linha;
  	int coluna;
  	int i;
  	int somaprod;

  	for(linha = 0; linha < a->lin; linha++){
    	for(coluna = 0; coluna < b->col; coluna++){
      		somaprod=0;
      		for(i = 0; i < a->col; i++){
				somaprod+=a->m[linha][i]*b->m[i][coluna];
			}
			resultante->m[linha][coluna]=somaprod;
		}
	}
	return resultante;
}

Matrix * get_rotation_matrix(double theta){

	//matriz composta pelos senos e cossenos do angulo theta
	Matrix * resultante = create_identity(3);

	resultante->m[1][1] = cos(theta);
	resultante->m[1][2] = -sen(theta);
	resultante->m[2][1] = sen(theta);
	resultante->m[2][2] = cos(theta);
		
	//return create_identity(3);
	return resultante;
}

Matrix * get_scale_matrix(double k){

	//uma matriz identidade onde a diagonal principal é composta pelo valor escalar
	Matrix * resultante = create_identity(3);
	resultante->m[1][1] = k;
	resultante->m[2][2] = k;

	//return create_identity(3);
	return resultante;
}

//a translação deve ser sempre a ultima a ser aplicada!!!!!!!
Matrix * get_translation_matrix(Vector v){

	Matrix * resultante = create_identity(3);
	resultante->m[1][3] = v.x;
	resultante->m[2][3] = v.y;

	//return create_identity(3);
	return resultante;
}

Matrix * get_transformation_matrix(Vector e1, Vector e2, Vector t){

	// TODO: implementar!
	Matrix bases = create_identity(3); //o que fazer com os dois vetores

	Matrix translacao = get_translation_matrix(t);
	Matrix * resultante = multiply(bases ,translacao);
		
	return create_identity(3);
}

Vector transform(Matrix * m, Vector v){

	Matrix * vetor = create_matrix(3, 1);

	vetor->m[1][1] = v.x;
	vetor->m[2][1] = v.y;
	vetor->m[3][1] = 1;

	Matrix * resultante = multiply(m, vetor);

	Vector * retorno =  create_vector();
	retorno.x = resultante->m[1][1];
	retorno.y = resultante->m[2][1];

	//return v;
	return retorno;
}

Matrix * get_observer_matrix(Vector position, Vector direction){

	// TODO: implementar!
		
	return create_identity(3);
}

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Programa principal. Le um arquivo de entrada com a seguinte estrutura: //
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

int main(int argc, char ** argv){

	////////////////
	//            //
	// Variaveis: //
	//            //
	////////////////

	char command[32];	
	char * input_file_name;		// nome do arquivo de entrada (com as definicoes do cena/desenho)
	char * output_file_name;	// nome do arquivo de saida (imagem gerada a partir das definicoes do arquivo de entrada)
	
	int width;			// largura da imagem a ser gerada
	int height;			// altura da imagem a ser gerada
	int background_color;		// cor de fundo da imagem a ser gerada
	
	int n_shapes;			// quantidade de shapes definidos no arquivo de entrada
	int i, j;			// variaveis de uso geral em iteracoes 
	
	Vector observer;		// vetor que representa a posicao do observador na cena a ser desenhada
	Vector direction;		// vetor que indica a direcao para a qual o observador olha na cena
	Vector v;			// variavel do tipo vetor de uso geral
	
	Shape * shapes;			// vetor de shapes. Usado para armazenar os shapes definidos no arquivo de entrada

	FILE * in;			// FILE handler para o arquivo de entrada
	Image * img;			// imagem na qual as operacoes de desenho serao realizada

	///////////////////////////////////////////
	//                                       //
	// Programa principal propriamente dito: //
	//                                       //
	///////////////////////////////////////////
	

	// verificacao dos parametros da linha de comando:
	
	if(argc != 3){

		printf("uso: %s <input_file_name> <output_file_name>\n", argv[0]);
		return 1;
	}
		
	input_file_name = argv[1]; 
	output_file_name = argv[2];
	
	// abertura do arquivo de entrada, e leitura dos parametros fixos (parametros da imagem e do observador, quantidade de shapes):

	in = fopen(input_file_name, "r");
	fscanf(in, "%d %d %d", &width, &height, &background_color);
	fscanf(in, "%lf %lf %lf %lf", &observer.x, &observer.y, &direction.x, &direction.y);
	fscanf(in, "%d", &n_shapes);
	
	img = create_imagem(width, height, background_color);	// criacao da imagem
	shapes = (Shape *) malloc (n_shapes * sizeof(Shape));	// alocacao do vetor de shapes com o tamanho adequado
	
	// leitura dos shapes definidos no arquivo de entrada:

	for(i = 0; i < n_shapes; i++){

		fscanf(in, "%d", &shapes[i].n);
		
		for(j = 0; j < shapes[i].n; j++) { 

			fscanf(in, "%lf %lf", &v.x, &v.y);
			shapes[i].vertices[j] = v;
		}
	}

	// leitura dos comandos de desenho, até que uma linha com o comando "END" seja encontrada:

	while(fscanf(in, "%s", command) == 1 && strcmp(command, END) != 0){

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

		int shape_id;		
		int color;		
		Vector t;
		Shape s;

		if(strcmp(command, DRAW_SHAPE) == 0 || strcmp(command, DRAW_SHAPE_BASE) == 0){
			
			if(strcmp(command, DRAW_SHAPE) == 0){

				double rotation;
				double scale;

				fscanf(in, "%d %d %lf %lf %lf %lf", &shape_id, &color, &rotation, &scale, &t.x, &t.y);

				// TODO: fala algo por aqui...
			}
			
			if(strcmp(command, DRAW_SHAPE_BASE) == 0){

				Vector e1, e2;
			
				fscanf(in, "%d %d %lf %lf %lf %lf %lf %lf", &shape_id, &color, &e1.x, &e1.y, &e2.x, &e2.y, &t.x, &t.y);

				// TODO: fala algo por aqui...
			}

			s = shapes[shape_id];

			for(i = 0; i < s.n - 1; i++){

				// TODO: fala algo por aqui...

				Vector v1 = s.vertices[i];
				Vector v2 = s.vertices[i + 1];
				draw_line(img, v1, v2, color);
			}
		}
		else {
			printf("Unknown command: '%s'\n", command);
		}
	}

	save_image(img, output_file_name);	// salva imagem no arquivo de saida
	free_image(img);			// libera os recursos de memoria usados pela imagem
	
	return 0;	
}
