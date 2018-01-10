from flask import Flask, send_from_directory, render_template
from flask import request
from flask import url_for
from datetime import datetime
from pytz import timezone


app = Flask(__name__)
@app.route("/")
def indexpage():
	now = datetime.now(timezone('Europe/London'))
	return render_template("index.html", time = now)
	#return "The time is {}".format(now)

@app.route("/family_table")
def family_table():
    return render_template("family_table.html")	
	
@app.route("/distribution")
def distribution():
    return "dis"	
	
@app.route("/AA_seq_list")
def AA_seq_list():
	newlist = 600+700
	return str(newlist)
	
@app.route("/relationship_AA")
def relationship_AA():
    return "rel"	
	
@app.route("/peptide_seq_ident")
def peptide_seq_ident():
    return "pep seq"	

@app.route("/upload_peptide")
def upload_peptide():
    return "upload"	
	
@app.route("/time2")
def time2():
    
    return "The time is2 {}".format(now)

@app.route("/profile/<name>")
def profile(name):
    return render_template("profile2.html", name=name)


@app.route('/user/<username>')
def show_user_profile(username):
	return render_template("style.html", username=username)
    # show the user profile for that user
    #return 'User %s' % username
	
@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        do_the_login()
    else:
        show_the_login_form()

@app.route("/visualise")
def visualise():
    return "lolnope"
	

	
#@app.route("/static")	
#def static1():
#	url_for('static', filename='style.css')
	
@app.route('/<path:path>')
def static_file(path):
    return app.send_static_file(path)

#if __name__ == "__main__":
#	app.run(debug=True)
