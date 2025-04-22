%%capture
!pip install nmrglue
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import sys

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from google.colab import drive
drive.mount('/content/drive')
from google.colab import drive
sys.path.append('/content/drive/MyDrive/py')
import ssnmr_proc.bruker as br
from scipy.io import wavfile
import scipy.signal as ss
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from google.colab import drive
import ssnmr_proc.bruker as br

# Adicione o caminho aos seus dados
# Adicione o caminho aos seus dados
sys.path.append('/content/drive/MyDrive/Mestrado/Resultados de pesquisa/RMN/dia-Cd(Im)2_1.3mm/9/pdata/1')
data_dir = r'/content/drive/MyDrive/Mestrado/Resultados de pesquisa/RMN/dia-Cd(Im)2_1.3mm/9/pdata/1'

# Carregue os dados
spc = br.ProcData(data_dir)

# Normaliza os dados
Z = spc.data / spc.data.real.max()
y = spc.ppm_scale
x = spc.ppm_scale_1

# Verifica a dimensão de x
print(x.ndim)

# Definir o índice da linha desejada
line_index = 529

# Cria a figura e os subplots
fig, axs = plt.subplots(2, 2, gridspec_kw={'height_ratios': [1, 5], 'width_ratios': [5, 1]}, figsize=(10, 8))

# Parâmetros do contorno
contour_start = 0.04           # valor inicial do nível de contorno
contour_num = 10              # número de níveis de contorno
contour_max = 1
contour_factor = np.exp(np.log(contour_max / contour_start) / contour_num)  # fator de escala entre os níveis de contorno

# Calcula os níveis de contorno
cl = contour_start * contour_factor ** np.arange(contour_num)

# Projeção dos dados ao longo do eixo 0 para o gráfico superior
proj_x = Z.sum(axis=0)
proj_x = proj_x / proj_x.max()

# Projeção dos dados ao longo do eixo 1 para o gráfico à direita
proj_y = Z.sum(axis=1)
proj_y = proj_y / proj_y.max()

# Gráfico superior (projeção x)
axs[0, 0].plot(x, proj_x, linestyle='solid', linewidth=1, color='k')
axs[0, 0].set_xlim([12, -3.5])
axs[0, 0].spines['top'].set_visible(False)
axs[0, 0].spines['right'].set_visible(False)
axs[0, 0].spines['bottom'].set_visible(False)
axs[0, 0].spines['left'].set_visible(False)
axs[0, 0].xaxis.set_visible(False)
axs[0, 0].yaxis.set_visible(False)
axs[0, 0].set_ylim([0, 1])

# Gráfico de contorno
contour = axs[1, 0].contour(x, y, Z, cl, cmap='winter')
axs[1, 0].set_ylim([19, -2])
axs[1, 0].set_xlim([12, -3.5])
axs[1, 0].set_xlabel(r'$^1$H-SQ /ppm', fontsize=16)
axs[1, 0].set_ylabel(r'$^1$H-DQ /ppm', fontsize=16)
axs[1, 0].plot(x, 2 * x, linestyle='dashed', linewidth=1.5, color='r')

# Gráfico superior direito (vazio)
axs[0, 1].spines['top'].set_visible(False)
axs[0, 1].spines['right'].set_visible(False)
axs[0, 1].spines['bottom'].set_visible(False)
axs[0, 1].spines['left'].set_visible(False)
axs[0, 1].xaxis.set_visible(False)
axs[0, 1].yaxis.set_visible(False)

# Gráfico à direita (projeção y)
axs[1, 1].plot(proj_y, y, linestyle='solid', linewidth=1, color='k')
axs[1, 1].set_ylim([19, -2])
axs[1, 1].spines['top'].set_visible(False)
axs[1, 1].spines['right'].set_visible(False)
axs[1, 1].spines['bottom'].set_visible(False)
axs[1, 1].spines['left'].set_visible(False)
axs[1, 1].xaxis.set_visible(False)
axs[1, 1].yaxis.set_visible(False)
axs[1, 1].set_xlim([0, 1])

# Adicionar eixo adicional para sobrepor o gráfico do slice
ax2 = axs[1, 0].twinx()
slice_data = Z[line_index, :]  # Dados do slice
scaled_slice_data = 18.5 + slice_data * (-20)  # Ajuste a escala conforme necessário
ax2.plot(x, scaled_slice_data, linestyle='solid', linewidth=2, color='purple')
ax2.set_ylim([19, -2])
ax2.set_yticks([])  # Remove as marcações do eixo y

# Ajuste dos espaçamentos
plt.subplots_adjust(hspace=0.05, wspace=0.05)

# Salva o gráfico
plt.savefig('grafico_comparacao.png', dpi=600)
plt.show()
