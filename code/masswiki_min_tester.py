# masswiki_min_tester.py

import requests

BASE = "https://masswiki.us-west-2.elasticbeanstalk.com/analysis/get_data"
ACCESS_TOKEN = "eyJraWQiOiJoeCtPbm1BUmpUWG0rZnNzZnptYVR2b2RIeG1Ra0dKbGVzc1hsZG5oTG5nPSIsImFsZyI6IlJTMjU2In0.eyJzdWIiOiI4ODYxYzM2MC1lMGExLTcwOTMtM2JjNC0wZDkxZDlmOGFkN2UiLCJjb2duaXRvOmdyb3VwcyI6WyJtYXNzd2lraS1sYWIiLCJtYXNzd2lraS1jb21tdW5pdHkiLCJVc2VycyJdLCJlbWFpbF92ZXJpZmllZCI6dHJ1ZSwiaXNzIjoiaHR0cHM6XC9cL2NvZ25pdG8taWRwLnVzLXdlc3QtMi5hbWF6b25hd3MuY29tXC91cy13ZXN0LTJfR2p0Y00wUENwIiwiY29nbml0bzp1c2VybmFtZSI6Ijg4NjFjMzYwLWUwYTEtNzA5My0zYmM0LTBkOTFkOWY4YWQ3ZSIsIm9yaWdpbl9qdGkiOiI3ZWFiN2VkOS1mMjA2LTRmMGMtODc3OC1kOWE4MjQ4MWQxMTUiLCJhdWQiOiIzaGdvczNhdGQxZWwxNmx0aDYxN2lpY29hbCIsImV2ZW50X2lkIjoiZThjMDdlM2ItMzc4Yi00NmY0LWEwMGItYTQ4NDgyNTk1YWE5IiwidG9rZW5fdXNlIjoiaWQiLCJhdXRoX3RpbWUiOjE3NTc2MzI4NDEsIm5hbWUiOiJaaXl1ZSBZYW5nIiwiZXhwIjoxNzU3NjM2NDQxLCJpYXQiOjE3NTc2MzI4NDEsImZhbWlseV9uYW1lIjoiWWFuZyIsImp0aSI6IjYwNmRkZTQzLTgwZTItNGFlNi04MTYzLWYwYzJhMmFmNTk4OCIsImVtYWlsIjoienl6eWFuZ0B1Y2RhdmlzLmVkdSJ9.LvlXkaEBslBJHCgfAftonrEfh30BeyZUDk48omMyIAfWQlQ8om5fL9W_0Nv80Kv_czf7niA7p29JNYhi0BA3NDk4aVrpVO09DyEnNqcHsJ2J7d242iB3Rp2gBoQe-XsQ_KoGMd_2zVw6SHE690kCSzvPbJIsUYPFvMYfdYrKaQ3eHTV_azxtgf-afdGqqVgA_jXpUjC24f9PHDsXNOyDkTqAoYylT7UhjZDejyry7_YqFYBt_0QFeyWquxQI40t8OLATw1yTN85EyQmUoRGGd9tV_Wu5KLYUNeXThXpEIvPqc-Dtb5wPcRbFAKRLkPMymCDh5AmWq9t3ZBu2JijP1A"  # use your corrected token here

def fetch_one(wiki_id: str):
    params = {
        "wiki_id": wiki_id,      # raw; DO NOT pre-encode
        "source": "binbase",
        "isPublic": "false",
    }
    headers = {
        "Accept": "application/json",
        "Authorization": f"Bearer {ACCESS_TOKEN}",
    }
    r = requests.get(BASE, params=params, headers=headers, timeout=20)
    print("URL sent:", r.url)      # confirm %2F (not %252F)
    print("Status:", r.status_code)
    print("Body (first 500):", r.text[:500])
    r.raise_for_status()
    return r.json()

# Try a couple of IDs
if __name__ == "__main__":
    fetch_one("aSW2RH3/5086")
    fetch_one("aSW2RH3/21")
