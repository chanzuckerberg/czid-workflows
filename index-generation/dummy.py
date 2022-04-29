import sys
import boto3
import csv
import mysql.connector

client = boto3.client('ssm')


def get_param(name: str) -> str:
    response = client.get_parameter(
        Name=name,
        WithDecryption=True,
    )
    return response["Parameter"]["Value"]


ENV = sys.argv[1]

mydb = mysql.connector.connect(
  host=get_param(f"/idseq-{ENV}-web/RDS_ADDRESS"),
  user=get_param(f"/idseq-{ENV}-web/DB_USERNAME"),
  password=get_param(f"/idseq-{ENV}-web/db_password"),
  port=int(get_param(f"/idseq-{ENV}-web/DB_PORT")),
  database=f"idseq_{ENV}",
)

mycursor = mydb.cursor()

mycursor.execute("CREATE TABLE taxon_lineages_new LIKE taxon_lineages")

myresult = mycursor.fetchall()

for x in myresult:
    print(x, file=sys.stderr)
