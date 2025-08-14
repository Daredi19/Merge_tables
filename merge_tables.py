# IMPORTS
import polars as pl
import os
import re
from rich import print
from rich.console import Console
from pathlib import Path
import csv
from itertools import islice
import argparse

console = Console()

# FUNCIONES

def get_df(file_path):
  '''Carga un archivo, detecta su separador y devuelve un dataframe.

  Args:
      file_path(str): Ruta del archivo.

  Returns:
      pl.DataFrame: DataFrame con los datos del archivo.
  '''

  try:

    file_path = Path(file_path).resolve()
    print(f"[yellow]:warning:[/] Se utilizará el siguiente archivo [bold green]{file_path.name}[/]")

    with open(file_path, 'r', newline='', encoding='utf-8') as f:

        # Lee (solo) las 5 primeras lineas y las concatena en una cadena ('' hace que no se meta nada en el medio) (Problema csvfile.readlines(): leería todo el archivo y despues se cogerian las 5 primeras)
        sample = ''.join(islice(f, 5)) 

        # La clase sniffer detecta el formato de un csv. Sniff es un metodo que devuelve el objeto Dialect. delimiter es un atributo de dialect que contiene el separador
        delimiter = csv.Sniffer().sniff(sample).delimiter 
        print(f"[yellow]:warning:[/] Se ha detectado el separador: [bold yellow]{repr(delimiter)}[/]")

        # Para devolver el puntero al principio del archivo
        f.seek(0) 
        
        # Polars lee todo el archivo y devuelve el dataframe, aumentamos el infer_schema_length para que detecte bien los tipos de dato de la columna prob
        df = pl.read_csv(f, separator=delimiter, infer_schema_length=10000) 

    return df

  except FileNotFoundError:
      console.print(f":x: Archivo no encontrado en la ruta: {file_path}")



def validate_annotations_table(df):
    '''Comprueba si la tabla de anotaciones que se proporciona tiene 15 o 35 columnas.
    Para ambos casos valida las tablas para comprobar que tienen un formato correcto.
    En función de cada caso, renombra y crea las columnas indicadas.

    Args:
        df(pl.DataFrame): Tabla obtenida del archivo de anotaciones.

    Returns:
        pl.DataFrame: DataFrame validado y con las columnas formateadas en función del número de columnas iniciales.
    '''

    if df.width == 15:
        header = ["X ID", "GENENAME", "DESCRIPTION"]
        missing_columns = []

        console.print(f"El DataFrame {"df"} contiene [bold green]15[/] columnas.")

        # Comprueba si todas las columnas indicadas en el encabezado están presentes. Si faltan las añade a la lista missing_columns
        for column in header:
            if column not in df.columns:
                missing_columns.append(column)    

        if len(missing_columns)!= 0:
            console.print(f"[bold red]:x:[/] La(s) columna(s) [italic purple]{missing_columns}[/] no está(n) presente(s) en el DataFrame {"df"}.")

        else:        
            console.print(f"[bold yellow]:warning: Se seleccionarán las siguientes columnas del DataFrame: [italic green]{header}[/]")
            console.print(f"[bold yellow]:warning: La columna [italic green] X ID[/], se modificará a [italic green] ID [/]")

            # Devuelve el df formateado con una nueva columna "ID" que contiene los tres primeros elementos que estaban separados por "::" en "X ID"            
            return df.select([pl.concat_str([pl.col("X ID").str.split("::").list.get(i) for i in range(3)], separator = "::").alias("ID"), 
                              pl.col("X ID").alias("old_ID"), 
                              "GENENAME",
                                "DESCRIPTION"])

    elif df.width == 35:
        header = ["ORF ID", "Gene name", "Gene length", "ORF length", "ORF start", "ORF end", 
                  "Strand", "Protein sequence", "Pfam", "InterPro", "GENENAME", "DESCRIPTION"]
        missing_columns = []

        console.print(f"El DataFrame de anotaciones contiene [bold green]35[/] columnas.")

        # Comprueba si todas las columnas indicadas en el encabezado están presentes. Si faltan las añade a la lista missing_columns
        for column in header:
            if column not in df.columns:
                missing_columns.append(column)    

        if len(missing_columns)!= 0:
            console.print(f"[bold red]:x:[/] La(s) columna(s) [italic purple]{missing_columns}[/] no está(n) presente(s) en el DataFrame {"df"}.")

        else:        
            console.print(f"[bold yellow]:warning: Se seleccionarán las siguientes columnas del DataFrame: [italic green]{header}[/]")
            return df.select([pl.col("ORF ID").alias("ID"), "Gene name", "Gene length", "ORF length", "ORF start", "ORF end", 
                  "Strand", "Protein sequence", "Pfam", "InterPro", "GENENAME", "DESCRIPTION"])
    
    else: 
        console.print(f"[bold red]:x:[/] DataFrame no válido. Debe contener 15 o 35 columnas.")



def validate_counts_table(df):
    '''Comprueba que la tabla de expresión diferencial (conteos) contiene las 6 columnas requeridas.
    Si las contiene renombra las columnas y crea una nueva con los ID simplificados..

    Args:
        df(pl.DataFrame): Tabla obtenida del archivo de expresión diferencial.

    Returns:
        pl.DataFrame: DataFrame validado y con las columnas formateadas.
    '''
    if df.width == 6:
        # Si el encabezado de la primera columna está vacío se le asigna el nombre de columna "old_ID"
        if df.columns[0] == '':
            df = df.rename({ '': "old_ID" }) #el fallo puede ser por esto?

        console.print(f"[bold yellow]:warning:[/] Se han cambiado los nombres de las columnas a: [italic green] old_ID, countsfiltered_ControlDMSO_mean, countsfiltered_DEHP_mean, theta, prob y log2FC[/].")
        console.print(f"[bold yellow]:warning:[/] Se ha generado una nueva columna:[italic green] ID [/].")

        # Crea un df nuevo con las columnas renombradas y la nueva columna "ID" que contiene la segunda parte del ID de "old_ID"
        df = df.select([
            pl.col(df.columns[0]).alias("old_ID"),
            pl.col("old_ID").str.slice(offset = pl.col("old_ID").str.find("~~") + 2).alias("ID"),
            pl.col(df.columns[1]).alias("countsfiltered_ControlDMSO_mean"),
            pl.col(df.columns[2]).alias("countsfiltered_DEHP_mean"),
            pl.col(df.columns[3]).alias("theta"),
            pl.col(df.columns[4]).alias("prob"),
            pl.col(df.columns[5]).alias("log2FC"),
            ])
        
        return df
    
    else:
        console.print(f"[bold red]:x:[/] DataFrame no válido. Debe contener 6 columnas.")
    


def merge_tables (annotations_path, counts_path, output_path, output_name):
    '''Une los dataframes de anotaciones y expresión diferencial, mediante la columna "ID".
    Se mantienen todas las filas de la tabla de expresión diferencial.
    Se guarda el resultado en un archivo csv
    
    Args:
        annotations_path(str): Ruta del archivo con la tabla de anotaciones.
        counts_path(str): Ruta del archivo con la tabla de expresión diferencial(conteos).
        output_path(str): Ruta donde se guardará el archivo csv generado.
        output_name(str): Nombre que recibe el archivo csv generado(especificar extensión).
    
    Returns:
        pl.DataFrame: Dataframe final con las tablas de anotaciones y expresión diferencial unidos.
    '''
    # Definir la ruta para el archivo de salida
    output_file = (f"{output_path}/{output_name}")

    # Cargar y validar la tabla de anotaciones
    annotations_df = validate_annotations_table(get_df(annotations_path))

    #Cargar y validar la tabla de conteos
    counts_df = validate_counts_table(get_df(counts_path))

    # Unir las tablas por la columna "ID" conservando todas las filas de la tabla de expresión diferencial
    final_df = counts_df.join(annotations_df, on = "ID", how = "left" )

    # Guarda el dataframe final en un archivo csv
    final_df.write_csv(output_file)  

    return final_df



def main(): 
    parser = argparse.ArgumentParser(description = "Une las tablas de anotaciones y de expresión diferencial")
    parser.add_argument("-a", "--annotations", required = True, help = "Ruta de la tabla de anotaciones")
    parser.add_argument("-c", "--counts", required = True, help = "Ruta a la tabla de expresión diferencial")
    parser.add_argument("-o", "--outputdir", required = True, help ="Ruta del directorio de salida")
    parser.add_argument("-n", "--name", required = True, help = "Nombre para el archivo de salida")

    args = parser.parse_args()

    merge_tables(
        annotations_path=args.annotations,
        counts_path=args.counts,
        output_path=args.outputdir,
        output_name=args.name
        )
    


# MAIN
if __name__ == "__main__":
    main()