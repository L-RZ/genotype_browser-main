from flask import Flask, jsonify, request, abort, render_template, redirect, url_for
from flask_compress import Compress
import imp, logging
import os
import re

from flask_login import LoginManager, login_user, UserMixin, login_required, logout_user
from werkzeug.security import generate_password_hash, check_password_hash
from wtforms import StringField, PasswordField
from wtforms.validators import DataRequired, EqualTo
from flask_wtf import FlaskForm

from utils import parse_chr, parse_region, ParseException, NotFoundException
from data import Datafetch
from search import Search

app = Flask(__name__, template_folder='../templates', static_folder='../static')
Compress(app)

app.secret_key = os.urandom(42)

login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'


#  TODO change the name and password before deploying
USERS = [
    {
        "id": 1,
        "name": 'lily',
        "password": generate_password_hash('123')
    },
    {
        "id": 2,
        "name": 'tom',
        "password": generate_password_hash('123')
    }
]

def get_user(user_name):
    for user in USERS:
        if user.get("name") == user_name:
            return user
    return None

class User(UserMixin):
    def __init__(self, user):
        self.username = user.get("name")
        self.password_hash = user.get("password")
        self.id = user.get("id")

    def verify_password(self, password):
        if self.password_hash is None:
            return False
        return check_password_hash(self.password_hash, password)

    def get_id(self):
        return self.id

    @staticmethod
    def get(user_id):
        if not user_id:
            return None
        for user in USERS:
            if user.get('id') == user_id:
                return User(user)
        return None

@login_manager.user_loader
def load_user(user_id):
    return User.get(user_id)

class LoginForm(FlaskForm):
    username = StringField('username', validators=[DataRequired()])
    password = PasswordField('password', validators=[DataRequired()])


@app.route('/login/', methods=('GET', 'POST'))
def login():
    form = LoginForm()
    emsg = None
    if form.validate_on_submit():
        user_name = form.username.data
        password = form.password.data
        user_info = get_user(user_name)
        if user_info is None:
            emsg = "username or password error"
        else:
            user = User(user_info)
            if user.verify_password(password):
                login_user(user, remember=True)
                return redirect(request.args.get('next') or url_for('index'))
            else:
                emsg = "username or password error"
    return render_template('login.html', form=form, emsg=emsg)


config = {}
try:
    _conf_module = imp.load_source('config', 'config.py')
except Exception as e:
    print('Could not load config.py')
    raise
config = {key: getattr(_conf_module, key) for key in dir(_conf_module) if not key.startswith('_')}

gunicorn_logger = logging.getLogger('gunicorn.error')
app.logger.handlers = gunicorn_logger.handlers
app.logger.setLevel(config['log_level'])

fetch = Datafetch(config)
search = Search(config)

@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
@login_required
def index(path):
    # return 'Hello World '
    return render_template('index.html')


@app.route('/logout')   # logout
@login_required
def logout():
    logout_user()
    return redirect(url_for('login'))


@app.route('/api/v1/find/<query>')
@login_required
def find(query):
    try:
        data_type = request.args.get('data_type')
        result = search.search(query, data_type)
    except ParseException as e:
        abort(400, 'could not parse given query to anything useful')
    except NotFoundException as e:
        abort(404, 'not found')
    return jsonify(result)


@app.route('/api/v1/variants/<variants>')
@login_required
def variants(variants):
    #print(request.args.to_dict())
    try:
        data_type = request.args.get('data_type')
        data = fetch.get_variants(variants, request.args.to_dict(), data_type)
    except ParseException as e:
        abort(400, 'could not parse given variant(s)')
    except NotFoundException as e:
        abort(404, 'variant(s) not in data')
    return jsonify(data)

@app.route('/api/v1/gene_variants/<gene>')
@login_required
def gene_variants(gene):
    try:
        data_type = request.args.get('data_type')
        data = fetch.get_gene_variants(gene, data_type)
    except ParseException as e:
        abort(400, 'could not parse given gene')
    except NotFoundException as e:
        abort(404, 'gene not in data')
    return jsonify(data)


@app.route('/api/v1/write_variants/<variants>')
@login_required
def write_variants(variants):
    try:
        data_type = request.args.get('data_type')
        status = fetch.write_variants(variants, request.args.to_dict(), data_type)
    except ParseException as e:
        abort(400, 'could not parse given variant(s)')
    except NotFoundException as e:
        abort(404, 'variant(s) not in data')
    return status


@app.route('/api/v1/range/<range>')
@login_required
def range(range):
    data_type = request.args.get('data_type')
    try:
        chr, start, end = parse_region(range)
        data = fetch.get_genomic_range_variants(chr, start, end, data_type)
    except ParseException as e:
        abort(400, 'could not parse given genomic range')
    except NotFoundException as e:
        abort(404, 'genomic range not in data')
    return jsonify(data)

# @app.route('/api/v1/clusterplot/<plot_type>/<variant>')
# def clusterplot(plot_type, variant):
#     try:
#         var = re.sub('-', '_', variant)
#         arr = var.split('_')
#         arr[0] = 'X' if arr[0] == '23' else arr[0]
#         exists_in_chip = fetch.check_var_in_chip(variant)
#         filename = config['cluster_plots_location'] + '/' + plot_type + '/' + '_'.join(arr) + '.png'
#         with open(filename, 'rb') as f:
#             blob = f.read()
#     except ParseException as e:
#         abort(400, 'could not parse given variant')
#     except FileNotFoundError as e:
#         if exists_in_chip:
#             abort(404, 'varaint exists in raw chip but no plot was found')
#         else:
#             abort(410, 'varaint does not exist in raw chip and no plot was found')
#     return blob


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080)
    # app.run(debug=True, port=8080, host='0.0.0.0')
