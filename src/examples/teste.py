import ssl

certs = list(ssl.enum_certificates("ROOT"))

print(f"{len(certs)} certificados encontrados\n")

ctx = ssl.create_default_context()

for i, cert in enumerate(certs):
    try:
        ctx.load_verify_locations(cadata=cert[0])
    except Exception as e:
        print(f"\nERRO no certificado {i}")
        print(e)
        print("\nTipo:", cert[1])
        print("Usos:", cert[2])
        break