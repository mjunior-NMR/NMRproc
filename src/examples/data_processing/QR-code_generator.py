# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 09:47:52 2025

@author: mcabe
"""

import qrcode

# Dados que você quer codificar no QR code
data = "https://www.exemplo.com"

# Criação do QR code
qr = qrcode.QRCode(
    version=1,
    error_correction=qrcode.constants.ERROR_CORRECT_L,
    box_size=10,
    border=4,
)

qr.add_data(data)
qr.make(fit=True)

# Criação da imagem do QR code
img = qr.make_image(fill_color="black", back_color="white")

# Salvando a imagem
img.save("qrcode_exemplo.png")
