import sys
import boto3
import mysql.connector

client = boto3.client('ssm')


def get_param(name: str) -> str:
    response = client.get_parameter(
        Name=name,
        WithDecryption=True,
    )
    return response["Parameter"]["Value"]


ENV = "sandbox"

mydb = mysql.connector.connect(
  host=get_param(f"/idseq-{ENV}-web/RDS_ADDRESS"),
  user=get_param(f"/idseq-{ENV}-web/DB_USERNAME"),
  password=get_param(f"/idseq-{ENV}-web/db_password"),
  port=int(get_param(f"/idseq-{ENV}-web/DB_PORT")),
  database=f"idseq_{ENV}",
)

mycursor = mydb.cursor()

mycursor.execute("SELECT * FROM taxon_lineages LIMIT 20")

myresult = mycursor.fetchall()

for x in myresult:
    print(x, file=sys.stderr)
